import groovy.json.JsonSlurper

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? file( params.multiqc_config, checkIfExists: true ) : []
ch_multiqc_logo            = params.multiqc_logo   ? file( params.multiqc_logo, checkIfExists: true ) : file("$projectDir/assets/mgnify_logo.png", checkIfExists: true)
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { SHORT_READS_QC                 } from '../subworkflows/local/short_reads_qc'
include { SHORT_READS_ASSEMBLY_QC        } from '../subworkflows/local/short_reads_assembly_qc'
include { SHORT_READS_ASSEMBLY_COVERAGE  } from '../subworkflows/local/short_reads_assembly_coverage'

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
include { SPADES                       } from '../modules/nf-core/spades/main'
include { MEGAHIT                      } from '../modules/nf-core/megahit/main'
include { QUAST                        } from '../modules/nf-core/quast/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SHORT_READS_ASSEMBLER {

    take:
    reads // tuple(meta), path(reads)

    main:

    ch_versions = Channel.empty()

    /***************************/
    /* Selecting the assembler */
    /***************************/
    /*
        The selection process ensures that:
        - The user selected assembler is always used (either from the samplesheet assembler column (with precedesnse) or the params.assembler)
        - Single-end reads are assembled with MEGAHIT, unless specified otherwise.
        - Paired-end reads are assembled with MetaSPAdes, unless specified otherwise
        - An error is raised if the assembler and read layout are incompatible (shouldn't happen...)
    */
    reads_by_assembler = reads.map { meta, reads ->
        def selected_assembler = meta.assembler;
        if ( selected_assembler == "megahit" || ( meta.single_end && selected_assembler == null ) ) {
            return [ meta + [assembler: "megahit", assembler_version: params.megahit_version], reads]
        } else if ( ["metaspades", "spades"].contains(selected_assembler) || ( !meta.single_end && selected_assembler == null ) ) {
            def xspades_assembler = selected_assembler ?: "metaspades" // Default to "metaspades" if the user didn't select one
            return [ meta + [assembler: xspades_assembler, assembler_version: params.spades_version], reads]
        } else {
            error "Incompatible assembler and/or reads layout. We can't assembly data that is. Reads - single end value: ${meta.single_end}."
        }
    }

    FASTQC_BEFORE (
        reads_by_assembler
    )

    ch_versions = ch_versions.mix(FASTQC_BEFORE.out.versions)

    SHORT_READS_QC(
        reads_by_assembler,
        params.reference_genome
    )

    FASTQC_AFTER (
        SHORT_READS_QC.out.qc_reads
    )

    /******************************************/
    /*  Reads that fail the following rules:  */
    /*  - Reads discarded by fastp > 90% (default value) */
    /*  - Less than 1k reads                  */
    /******************************************/
    extended_qc = SHORT_READS_QC.out.fastp_json.map { meta, json -> {
            json_txt = new JsonSlurper().parseText(json.text)
            bf_total_reads = json_txt?.summary?.before_filtering?.total_reads ?: 0;
            af_total_reads = json_txt?.summary?.after_filtering?.total_reads ?: 0;
            reads_qc_meta = [
                "low_reads_count": af_total_reads <= params.low_reads_count_threshold,
                "filter_ratio_threshold_exceeded": af_total_reads == 0 || ((af_total_reads / bf_total_reads) <= params.filter_ratio_threshold )
            ]
            return [meta, reads_qc_meta]
        }
    }

    extended_reads_qc = SHORT_READS_QC.out.qc_reads.join( extended_qc )

    extended_reads_qc.branch { meta, reads, reads_qc_meta ->
        // Filter out failed reads //
        qc_failed: reads_qc_meta.low_reads_count || reads_qc_meta.filter_ratio_threshold_exceeded
        megahit: meta.assembler == "megahit"
        xspades: ["metaspades", "spades"].contains(meta.assembler)
    }.set { qc_filtered_reads }

    ch_versions = ch_versions.mix(SHORT_READS_QC.out.versions)

    /*********************/
    /*     Assembly     */
    /********************/
    SPADES(
        qc_filtered_reads.xspades.map { meta, reads, _ -> [meta, reads, [], []] },
        [], // yml input parameters, which we don't use
        []  // hmm, not used
    )

    ch_versions = ch_versions.mix(SPADES.out.versions)

    MEGAHIT(
        qc_filtered_reads.megahit.map { meta, reads, _ -> [meta, reads] }
    )

    assembly = SPADES.out.contigs.mix( MEGAHIT.out.contigs )

    ch_versions = ch_versions.mix(MEGAHIT.out.versions)

    // Clean the assembly contigs //
    SHORT_READS_ASSEMBLY_QC(
        assembly,
        params.reference_genome
    )

    ch_versions = ch_versions.mix(SHORT_READS_ASSEMBLY_QC.out.versions)

    // Coverage //
    SHORT_READS_ASSEMBLY_COVERAGE(
        SHORT_READS_ASSEMBLY_QC.out.filtered_contigs.join( SHORT_READS_QC.out.qc_reads, remainder: false )
    )

    ch_versions = ch_versions.mix(SHORT_READS_ASSEMBLY_COVERAGE.out.versions)

    // Stats //
    /* The QUAST module was modified to run metaQUAST instead */
    QUAST(
        SHORT_READS_ASSEMBLY_QC.out.filtered_contigs,
        [ [], [] ], // reference
        [ [], [] ]  // gff
    )

    ch_versions = ch_versions.mix(QUAST.out.versions)

    emit:
    fastqc_before_zip                    = FASTQC_BEFORE.out.zip                                // tuple(meta)
    qc_failed                            = qc_filtered_reads.qc_failed                          // tuple(meta)
    fastqc_after_zip                     = FASTQC_AFTER.out.zip                                // tuple(meta)
    assembly_coverage_samtools_idxstats  = SHORT_READS_ASSEMBLY_COVERAGE.out.samtools_idxstats // tuple(meta)
    quast_results                        = QUAST.out.results                                   // tuple(meta)
    versions                             = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
