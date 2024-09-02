/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsSummaryLog; paramsSummaryMap; samplesheetToList } from 'plugin/nf-schema'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

validateParameters()

if (params.help) {
    log.info paramsHelp("nextflow run ebi-metagenomics/longreadsassembly --help")
    exit 0
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

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
// include { PACBIO_LQ       } from '../subworkflows/local/pacbio_lq'
// include { PACBIO_HIFI     } from '../subworkflows/local/pacbio_hifi'

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
def multiqc_report = []

workflow LONGREADSASSEMBLY {

    ch_versions = Channel.empty()
    longReads = Channel.empty()
    fetch_tool_metadata = Channel.empty()

    if ( params.samplesheet ) {

        longReads = { study_accession, reads_accession, fq1, library_layout, library_strategy, assembler, assembler_config, assembly_memory ->
            return tuple(
                [
                    "id": reads_accession,
                    "study_accession": study_accession,
                    "library_strategy": library_strategy,
                    "library_layout": library_layout,
                    "single_end": true,
                    "assembler": assembler ?: params.assembler,
                    "assembler_config": assembler_config ?: params.assembler_config,
                    "assembly_memory": assembly_memory ?: params.assembly_memory
                ],
                [fq1]
            )
        }

        samplesheet = Channel.fromList(samplesheetToList(params.samplesheet, "./assets/schema_input.json"))

        fetch_reads_transformed = samplesheet.map(longReads)

    } else {
        // TODO: remove when the fetch tools gets published on bioconda
        fetch_tool_config = file("${projectDir}/assets/fetch_tool_anonymous.json", checkIfExists: true)

        if ( params.private_study ) {
            fetch_tool_config = file("${projectDir}/assets/fetch_tool_credentials.json", checkIfExists: true)
        }

        FETCHTOOL_READS(
            [ [id: params.reads_accession], params.study_accession, params.reads_accession ],
            fetch_tool_config
        )

        ch_versions = ch_versions.mix(FETCHTOOL_READS.out.versions)

        // Push the library strategy into the meta of the reads, this is to make it easier to handle downstream
        fetch_reads_transformed = FETCHTOOL_READS.out.reads.map { meta, reads, library_strategy, library_layout, platform -> {
                [ meta + [
                    //  -- The metadata will be overriden by the parameters -- //
                    "assembler": params.assembler,
                    "assembly_memory": params.assembly_memory,
                    "assembler_config": params.assembler_config,
                    "library_strategy": params.library_strategy ?: library_strategy,
                    "library_layout": params.library_layout ?: library_layout,
                    "single_end": params.single_end ?: library_layout == "single",
                    "platform": params.platform ?: platform
                ], reads ]
            }
        }

        // Metadata for MultiQC
        fetch_tool_metadata = FETCHTOOL_READS.out.metadata_tsv.map { it[1] }.collectFile(
            name: 'fetch_tool_mqc.tsv',
            newLine: true,
            keepHeader: true,
            skip: 1
        )
    }

    LONG_READS_QC (
        fetch_reads_transformed, 
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
            if (params.assembler_config == "nano-raw" || meta.quality == "low") {
                return [meta + ["assembler_config": "nano-raw"], reads]
            } else if (params.assembler_config == "nano-hq" || meta.quality == "high") {
                return [meta + ["assembler_config": "nano-hq"], reads]
            }
        } else if (meta.platform == "pacbio") {
            if (params.assembler_config == "pacbio-raw" || meta.quality == "low") {
                return [meta + ["assembler_config": "pacbio-raw"], reads]
            } else if (params.assembler_config == "pacbio-hifi" || meta.quality == "high") {
                return [meta + ["assembler_config": "pacbio-hifi"], reads]
            }
        } else {
            error "Incompatible configuration"
        }
    }
    
    reads_assembler_config.branch { meta, reads ->
        lq_ont: meta.assembler_config == "nano-raw"
        hq_ont: meta.assembler_config == "pacbio-raw"
        lq_pacbio: meta.assembler_config == "nano-hq"
        hq_pacbio: meta.assembler_config == "pacbio-hifi"
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

    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )

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
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
