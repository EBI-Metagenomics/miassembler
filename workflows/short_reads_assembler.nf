/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { DOWNLOAD_FROM_FIRE            } from '../modules/local/download_from_fire.nf'

include { SHORT_READS_QC                } from '../subworkflows/local/short_reads_qc'
include { SHORT_READS_ASSEMBLY_QC       } from '../subworkflows/local/short_reads_assembly_qc'
include { SHORT_READS_ASSEMBLY_COVERAGE } from '../subworkflows/local/short_reads_assembly_coverage'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC as FASTQC_BEFORE       } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_AFTER        } from '../modules/nf-core/fastqc/main'
include { SPADES                        } from '../modules/nf-core/spades/main'
include { MEGAHIT                       } from '../modules/nf-core/megahit/main'
include { QUAST                         } from '../modules/nf-core/quast/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SHORT_READS_ASSEMBLER {
    take:
    input_reads // tuple(meta), path(reads)

    main:

    def ch_versions = Channel.empty()
    def reads_to_assemble = input_reads

    // If running for a private study on EBI infrastructure //
    if (params.private_study) {
        /*
         * For private studies we need to bypass Nextflow S3 integration until https://github.com/nextflow-io/nextflow/issues/4873 is fixed
         * The EBI parameter is needed as this only works on EBI network, FIRE is not accessible otherwise
        */
        DOWNLOAD_FROM_FIRE(
            input_reads
        )

        ch_versions = ch_versions.mix(DOWNLOAD_FROM_FIRE.out.versions.first())

        reads_to_assemble = DOWNLOAD_FROM_FIRE.out.reads
    }

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
    def reads_by_assembler = reads_to_assemble.map { meta, reads ->
        def selected_assembler = meta.assembler
        if (selected_assembler == "megahit" || (meta.single_end && selected_assembler == null)) {
            return [meta + [assembler: "megahit", assembler_version: params.megahit_version], reads]
        }
        else if (["metaspades", "spades"].contains(selected_assembler) || (!meta.single_end && selected_assembler == null)) {
            def xspades_assembler = selected_assembler ?: "metaspades"
            // Default to "metaspades" if the user didn't select one
            return [meta + [assembler: xspades_assembler, assembler_version: params.spades_version], reads]
        }
        else {
            error("Incompatible assembler and/or reads layout. We can't assembly data that is. Reads - single end value: ${meta.single_end}.")
        }
    }

    FASTQC_BEFORE(
        reads_by_assembler
    )
    ch_versions = ch_versions.mix(FASTQC_BEFORE.out.versions)

    SHORT_READS_QC(
        reads_by_assembler,
        params.reference_genome
    )
    ch_versions = ch_versions.mix(SHORT_READS_QC.out.versions)

    FASTQC_AFTER(
        SHORT_READS_QC.out.qc_reads
    )

    /******************************************/
    /*  Reads that fail the following rules:  */
    /*  - Reads discarded by fastp > 90% (default value) */
    /*  - Less than 1k reads                  */
    /******************************************/
    def extended_qc = SHORT_READS_QC.out.fastp_json.map { meta, json ->
        {
            def json_txt = new groovy.json.JsonSlurper().parseText(json.text)
            def bf_total_reads = json_txt.summary.before_filtering.total_reads ?: 0
            def af_total_reads = json_txt.summary.after_filtering.total_reads ?: 0
            def reads_qc_meta = [
                "low_reads_count": af_total_reads <= params.short_reads_low_reads_count_threshold,
                "filter_ratio_threshold_exceeded": af_total_reads == 0 || ((af_total_reads / bf_total_reads) <= params.short_reads_filter_ratio_threshold)
            ]
            return [meta, reads_qc_meta]
        }
    }

    def extended_reads_qc = SHORT_READS_QC.out.qc_reads.join(extended_qc)

    extended_reads_qc
        .branch { meta, _reads, reads_qc_meta ->
            qc_failed: reads_qc_meta.low_reads_count || reads_qc_meta.filter_ratio_threshold_exceeded
            megahit: meta.assembler == "megahit"
            xspades: ["metaspades", "spades"].contains(meta.assembler)
        }
        .set { qc_filtered_reads }

    /*********************/
    /*     Assembly     */
    /********************/
    SPADES(
        qc_filtered_reads.xspades.map { meta, reads, __ -> [meta, reads, [], []] },
        [],
        []
    )
    ch_versions = ch_versions.mix(SPADES.out.versions)

    MEGAHIT(
        qc_filtered_reads.megahit.map { meta, reads, __ -> [meta, reads] }
    )
    ch_versions = ch_versions.mix(MEGAHIT.out.versions)

    assembly = SPADES.out.contigs.mix(MEGAHIT.out.contigs)

    // Clean the assembly contigs //
    SHORT_READS_ASSEMBLY_QC(
        assembly,
        params.reference_genome
    )
    ch_versions = ch_versions.mix(SHORT_READS_ASSEMBLY_QC.out.versions)

    // Coverage //
    SHORT_READS_ASSEMBLY_COVERAGE(
        SHORT_READS_ASSEMBLY_QC.out.filtered_contigs.join(SHORT_READS_QC.out.qc_reads, remainder: false)
    )

    ch_versions = ch_versions.mix(SHORT_READS_ASSEMBLY_COVERAGE.out.versions)

    // Stats //
    /* The QUAST module was modified to run metaQUAST instead */
    QUAST(
        SHORT_READS_ASSEMBLY_QC.out.filtered_contigs,
        [[], []],
        [[], []]
    )

    ch_versions = ch_versions.mix(QUAST.out.versions)

    emit:
    fastqc_before_zip                   = FASTQC_BEFORE.out.zip // tuple(meta)
    qc_failed                           = qc_filtered_reads.qc_failed // tuple(meta)
    fastqc_after_zip                    = FASTQC_AFTER.out.zip // tuple(meta)
    assembly_coverage_samtools_idxstats = SHORT_READS_ASSEMBLY_COVERAGE.out.samtools_idxstats // tuple(meta)
    quast_results                       = QUAST.out.results // tuple(meta)
    versions                            = ch_versions
}
