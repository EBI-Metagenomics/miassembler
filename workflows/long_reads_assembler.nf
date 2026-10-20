/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { LONG_READS_QC                } from '../subworkflows/local/long_reads_qc'

include { ONT_LQ                       } from '../subworkflows/local/ont_lq'
include { ONT_HQ                       } from '../subworkflows/local/ont_hq'
include { PACBIO_LQ                    } from '../subworkflows/local/pacbio_lq'
include { PACBIO_HIFI                  } from '../subworkflows/local/pacbio_hifi'

include { LONG_READS_ASSEMBLY_QC       } from '../subworkflows/local/long_reads_assembly_qc'
include { FRAMESHIFT_CORRECTION        } from '../subworkflows/local/frameshift_correction'
include { LONG_READS_ASSEMBLY_COVERAGE } from '../subworkflows/local/long_reads_assembly_coverage'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC as FASTQC_BEFORE      } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_AFTER       } from '../modules/nf-core/fastqc/main'
include { QUAST                        } from '../modules/nf-core/quast/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary

workflow LONG_READS_ASSEMBLER {

    take:
    input_reads // tuple(meta), path(input_reads)

    main:

    ch_versions = channel.empty()

    // Long reads are always assembled with Flye
    def reads_flye = input_reads.map { meta, reads ->
        [meta + [assembler: "flye", assembler_version: params.flye_version], reads]
    }

    FASTQC_BEFORE (
        reads_flye
    )
    ch_versions = ch_versions.mix(FASTQC_BEFORE.out.versions)

    LONG_READS_QC (
        reads_flye
    )
    ch_versions = ch_versions.mix(LONG_READS_QC.out.versions)

    FASTQC_AFTER (
        LONG_READS_QC.out.qc_reads
    )
    ch_versions = ch_versions.mix(FASTQC_AFTER.out.versions)

    /*********************************************************************************/
    /* Selecting the combination of adapter trimming, assembler, and post-processing */
    /*********************************************************************************/
    /*
        The selection process ensures that:
        - The user selected assembler configuration is always used (either from the samplesheet assembler column (with precedence) or the params.assembler)
        - Low-quality ONT reads are trimmed with chopper and assembled with flye --nano-raw), unless specified otherwise.
        - High-quality ONT reads are trimmed with porechop_abi and assembled with flye --nano-hq), unless specified otherwise.
        - Low-quality pacbio reads are trimmed with chopper and assembled with flye --pacbio-raw), unless specified otherwise.
        - High-quality pacbio reads are trimmed with HiFiAdapterFilt and assembled with flye --pacbio-hifi), unless specified otherwise.
        Extra polishing steps are applied to low-quality reads. All subworkflows also apply post-assembly host decontamination.
    */

    // meta.assembler_config is never overriden if provided as input.
    // If no input was provided, quality and platform will determine the assembly type
    reads_assembler_config = LONG_READS_QC.out.qc_reads.map { meta, reads ->
        if (!meta.assembler_config) {
            if (meta.platform == "ont") {
                if (meta.quality == "low") {
                    return [meta + ["assembler_config": "nano-raw"], reads]
                } else if (meta.quality == "high") {
                    return [meta + ["assembler_config": "nano-hq"], reads]
                }
            } else if (meta.platform == "pb") {
                if (meta.quality == "low") {
                    return [meta + ["assembler_config": "pacbio-raw"], reads]
                } else if (meta.quality == "high") {
                    return [meta + ["assembler_config": "pacbio-hifi"], reads]
                }
            } else {
                error "Invalid quality ${meta.quality} or platform ${meta.platform} for ${meta.id}"
            }
        } else {
            return [meta, reads]
        }
    }

    reads_assembler_config.branch { meta, _reads ->
        lq_ont: meta.assembler_config == "nano-raw"
        hq_ont: meta.assembler_config == "nano-hq"
        lq_pacbio: meta.assembler_config == "pacbio-raw"
        hq_pacbio: meta.assembler_config == "pacbio-hifi"
        // For example nano-corr and pacbio-corr: not assembled, reported in qc_failed_runs.csv
        unsupported: true
    }.set {subworkflow_platform_reads}

    def unsupported_reads = subworkflow_platform_reads.unsupported.map { meta, reads ->
        [meta + ["unsupported_assembler_config": true], reads]
    }

    ONT_LQ(
        subworkflow_platform_reads.lq_ont
    )
    ch_versions = ch_versions.mix(ONT_LQ.out.versions)

    ONT_HQ(
        subworkflow_platform_reads.hq_ont
    )
    ch_versions = ch_versions.mix(ONT_HQ.out.versions)

    PACBIO_LQ(
        subworkflow_platform_reads.lq_pacbio
    )
    ch_versions = ch_versions.mix(PACBIO_LQ.out.versions)

    PACBIO_HIFI(
        subworkflow_platform_reads.hq_pacbio
    )
    ch_versions = ch_versions.mix(PACBIO_HIFI.out.versions)

    def assembly = ONT_LQ.out.contigs.mix(ONT_HQ.out.contigs,
                                        PACBIO_LQ.out.contigs,
                                        PACBIO_HIFI.out.contigs)

    // /**********************************************************************************/
    // /* Post-assembly: host decontamination, frame-shift correction, coverage and stats */
    // /**********************************************************************************/

    LONG_READS_ASSEMBLY_QC(assembly)
    ch_versions = ch_versions.mix(LONG_READS_ASSEMBLY_QC.out.versions)

    def decontaminated_assembly = LONG_READS_ASSEMBLY_QC.out.contigs

    decontaminated_assembly.branch { meta, _contigs ->
        lq: meta.quality == "low"
        hq: meta.quality == "high"
    }.set{low_high_quality_contigs}

    FRAMESHIFT_CORRECTION(
        low_high_quality_contigs.lq
    )
    ch_versions = ch_versions.mix(FRAMESHIFT_CORRECTION.out.versions)

    final_contigs = FRAMESHIFT_CORRECTION.out.corrected_contigs.mix(
                        low_high_quality_contigs.hq)

    def final_contigs_by_run = final_contigs.map { meta, contigs ->
        [meta.subMap("study_accession", "id", "platform", "assembler_config"), meta, contigs]
    }

    def reads_assembler_config_by_run = reads_assembler_config.map { meta, reads ->
        [meta.subMap("study_accession", "id", "platform", "assembler_config"), reads]
    }

    def final_contigs_reads = final_contigs_by_run
        .join(reads_assembler_config_by_run)
        .map { _run_meta, meta, contigs, reads ->
            [meta, contigs, reads]
        }

    LONG_READS_ASSEMBLY_COVERAGE(
        final_contigs_reads,
        LONG_READS_QC.out.fastp_json
    )
    ch_versions = ch_versions.mix(LONG_READS_ASSEMBLY_COVERAGE.out.versions)

    // Stats //
    /* The QUAST module was modified to run metaQUAST instead */
    QUAST(
        final_contigs,
        [ [], [] ], // reference
        [ [], [] ]  // gff
    )

    ch_versions = ch_versions.mix(QUAST.out.versions)

    emit:
    fastqc_before_zip                    = FASTQC_BEFORE.out.zip                              // tuple(meta)
    fastqc_after_zip                     = FASTQC_AFTER.out.zip                               // tuple(meta)
    qc_failed_all                        = unsupported_reads                                  // tuple(meta)
    assembly_coverage_samtools_idxstats  = LONG_READS_ASSEMBLY_COVERAGE.out.samtools_idxstats // tuple(meta)
    quast_results                        = QUAST.out.results                                  // tuple(meta)
    decontamination_mqc                  = LONG_READS_ASSEMBLY_QC.out.decontamination_mqc     // tuple(meta), path(mqc_tsv)
    decontaminated_contigs               = LONG_READS_ASSEMBLY_QC.out.contigs                 // tuple(meta), path(assembly_fasta)
    contaminated_contigs                 = LONG_READS_ASSEMBLY_QC.out.contaminated_contigs    // tuple(meta), val(step), path(tsv)
    versions                             = ch_versions

    // Published outputs //
    fastqc_before                        = FASTQC_BEFORE.out.html.mix(FASTQC_BEFORE.out.zip)
    fastqc_after                         = FASTQC_AFTER.out.html.mix(FASTQC_AFTER.out.zip)
    fastplong_reports                    = LONG_READS_QC.out.fastplong_reports
    contigs                              = final_contigs                                      // LQ (frame-shift corrected) + HQ
    decontamination                      = LONG_READS_ASSEMBLY_QC.out.decontamination
    coverage                             = LONG_READS_ASSEMBLY_COVERAGE.out.coverage_depth.mix(
                                               LONG_READS_ASSEMBLY_COVERAGE.out.assembly_coverage_json
                                           )
    quast                                = QUAST.out.results.mix(QUAST.out.tsv)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
