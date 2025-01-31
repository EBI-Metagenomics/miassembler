/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { FETCHTOOL_READS } from '../modules/local/fetchtool_reads'
include { LONG_READS_QC   } from '../subworkflows/local/long_reads_qc'

include { ONT_LQ          } from '../subworkflows/local/ont_lq'
include { ONT_HQ          } from '../subworkflows/local/ont_hq'
include { PACBIO_LQ       } from '../subworkflows/local/pacbio_lq'
include { PACBIO_HIFI     } from '../subworkflows/local/pacbio_hifi'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary

workflow LONG_READS_ASSEMBLER {

    take:
    reads // tuple(meta), path(reads)

    main:

    ch_versions = Channel.empty()

    LONG_READS_QC (
        reads,
        params.reference_genome
    )
    ch_versions = ch_versions.mix(LONG_READS_QC.out.versions)

    /*********************************************************************************/
    /* Selecting the combination of adapter trimming, assembler, and post-processing */
    /*********************************************************************************/
    /*
        The selection process ensures that:
        - The user selected assembler configuration is always used (either from the samplesheet assembler column (with precedence) or the params.assembler)
        - Low-quality ONT reads are trimmed with canu and assembled with flye --nano-corr/raw), unless specified otherwise.
        - High-quality ONT reads are trimmed with porechob_abi and assembled with flye --nano-hq), unless specified otherwise.
        - Low-quality pacbio reads are trimmed with canu and assembled with flye --pacbio-corr/raw), unless specified otherwise.
        - High-quality pacbio reads are trimmed with HiFiAdapterFilt and assembled with flye --pacbio-hifi), unless specified otherwise.
        Extra polishing steps are applied to low-quality reads. All subworkflows also apply post-assembly host decontamination.
    */

    reads_assembler_config = LONG_READS_QC.out.qc_reads.map { meta, reads ->
        if (meta.platform == "ont") {
            if (params.long_reads_assembler_config == "nano-raw" || meta.quality == "low") {
                return [meta + ["assembler_config": "nano-raw"], reads]
            } else if (params.long_reads_assembler_config == "nano-hq" || meta.quality == "high") {
                return [meta + ["assembler_config": "nano-hq"], reads]
            }
        } else if (meta.platform == "pacbio") {
            if (params.long_reads_assembler_config == "pacbio-raw" || meta.quality == "low") {
                return [meta + ["assembler_config": "pacbio-raw"], reads]
            } else if (params.long_reads_assembler_config == "pacbio-hifi" || meta.quality == "high") {
                return [meta + ["assembler_config": "pacbio-hifi"], reads]
            }
        } else {
            error "Incompatible configuration"
        }
    }

    /*********************************************************************************/
    /* Selecting the combination of adapter trimming, assembler, and post-processing */
    /*********************************************************************************/
    /*
        The selection process ensures that:
        - The user selected assembler configuration is always used (either from the samplesheet assembler column (with precedence) or the params.assembler)
        - Low-quality ONT reads are trimmed with canu and assembled with flye --nano-corr/raw), unless specified otherwise.
        - High-quality ONT reads are trimmed with porechob_abi and assembled with flye --nano-hq), unless specified otherwise.
        - Low-quality pacbio reads are trimmed with canu and assembled with flye --pacbio-corr/raw), unless specified otherwise.
        - High-quality pacbio reads are trimmed with HiFiAdapterFilt and assembled with flye --pacbio-hifi), unless specified otherwise.
        Extra polishing steps are applied to low-quality reads. All subworkflows also apply post-assembly host decontamination.
    */

    reads_assembler_config = LONG_READS_QC.out.qc_reads.map { meta, reads ->
        if (meta.platform == "ont") {
            if (params.long_reads_assembler_config == "nano-raw" || meta.quality == "low") {
                return [meta + ["long_reads_assembler_config": "nano-raw"], reads]
            } else if (params.long_reads_assembler_config == "nano-hq" || meta.quality == "high") {
                return [meta + ["long_reads_assembler_config": "nano-hq"], reads]
            }
        } else if (meta.platform == "pacbio") {
            if (params.long_reads_assembler_config == "pacbio-raw" || meta.quality == "low") {
                return [meta + ["long_reads_assembler_config": "pacbio-raw"], reads]
            } else if (params.long_reads_assembler_config == "pacbio-hifi" || meta.quality == "high") {
                return [meta + ["long_reads_assembler_config": "pacbio-hifi"], reads]
            }
        } else {
            error "Incompatible configuration"
        }
    }

    reads_assembler_config.branch { meta, reads ->
        lq_ont: meta.long_reads_assembler_config == "nano-raw"
        hq_ont: meta.long_reads_assembler_config == "pacbio-raw"
        lq_pacbio: meta.long_reads_assembler_config == "nano-hq"
        hq_pacbio: meta.long_reads_assembler_config == "pacbio-hifi"
    }.set {subworkflow_platform_reads}

    ONT_LQ(
        subworkflow_platform_reads.lq_ont
    )

    ONT_HQ(
        subworkflow_platform_reads.hq_ont
    )

    // PACBIO_LQ(
    //     subworkflow_platform_reads.lq_pacbio.map { meta, reads -> [meta, reads] }
    // )

    // PACBIO_HIFI(
    //     subworkflow_platform_reads.hq_pacbio.map { meta, reads -> [meta, reads] }
    // )

    assembly = ONT_LQ.out.contigs.mix( ONT_HQ.out.contigs )//, PACBIO_LQ.out.contigs, PACBIO_HIFI.out.contigs )

    /*************************************/
    /* Post-assembly: coverage and stats */
    /*************************************/

    //
    // MODULE: Run FastQC
    //
    // FASTQC (
    //     INPUT_CHECK.out.reads
    // )
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: MultiQC
    //
    // workflow_summary    = WorkflowLongreadsassembly.paramsSummaryMultiqc(workflow, summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    // methods_description    = WorkflowLongreadsassembly.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    // ch_methods_description = Channel.value(methods_description)

    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList()
    // )
    // multiqc_report = MULTIQC.out.report.toList()

    emit:
    versions                             = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
