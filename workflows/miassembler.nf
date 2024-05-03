/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

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
include { FETCHTOOL_READS    } from '../modules/local/fetchtool_reads'
include { FETCHTOOL_METADATA } from '../modules/local/fetchtool_metadata'
include { READS_QC           } from '../subworkflows/local/reads_qc'
include { ASSEMBLY_QC        } from '../subworkflows/local/assembly_qc'
include { ASSEMBLY_COVERAGE  } from '../subworkflows/local/assembly_coverage'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SPADES                      } from '../modules/nf-core/spades/main'
include { MEGAHIT                     } from '../modules/nf-core/megahit/main'
include { QUAST                       } from '../modules/nf-core/quast/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

def metaSorter = { a, b ->
    // Check if both a and b are LinkedHashMap
    if (a instanceof LinkedHashMap && b instanceof LinkedHashMap) {
        // Compare the lengths
        return b.size() <=> a.size()
    } else if (a instanceof LinkedHashMap) {
        // LinkedHashMaps (like meta) is considered bigger than value
        return 1
    } else if (b instanceof LinkedHashMap) {
        return -1
    } else {
        return a <=> b
    }
}

workflow MIASSEMBLER {

    ch_versions = Channel.empty()
    
    fetch_tool_config = file("$projectDir/assets/fetch_tool_anonymous.json")
    if ( params.private_study ) {
        fetch_tool_config = file("$projectDir/assets/fetch_tool_credentials.json")
    }

    // Download project metadata //
    FETCHTOOL_METADATA(
        [ [id: params.reads_accession], params.study_accession, params.reads_accession ],
        fetch_tool_config
    )
    
    ch_versions = ch_versions.mix(FETCHTOOL_METADATA.out.versions)
    
    // Download reads //
    FETCHTOOL_READS(
        [ [id: params.reads_accession], params.study_accession, params.reads_accession ],
        fetch_tool_config
    )

    ch_versions = ch_versions.mix(FETCHTOOL_READS.out.versions)

    FASTQC (
        FETCHTOOL_READS.out.reads
    )

    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // Perform QC on reads //
    READS_QC(
        FETCHTOOL_READS.out.reads,
        params.reference_genome
    )

    /*
    Single end reads // paired end reads distinction
        We need to split single-end and paired-end reads.
        Single-end reads are always assembled with MEGAHIT.
    */

    READS_QC.out.reads.branch { meta, reads ->
        xspades: ["metaspades", "spades"].contains(params.assembler)
            && meta.single_end == false
            || FETCHTOOL_METADATA.out.lib_strategy.contains("METATRANSCRIPTOMIC")
        megahit: params.assembler == "megahit" || meta.single_end == true
    }.set { qc_reads }

    ch_versions = ch_versions.mix(READS_QC.out.versions)

    /* Assembly */
    /* -- Clarification --
        At the moment, the pipeline only processes one set of reads at a time.
        Therefore, running Spades, metaSpades, or MEGAHIT are mutually exclusive.
        In order to support multiple runs, we need to refactor the code slightly.
        We will need to use `.join()` to keep assemblies together with their corresponding reads.
    */

    SPADES(
        qc_reads.xspades.map { meta, reads -> [meta, reads, [], []] },
        params.assembler,
        [], // yml input parameters, which we don't use
        []  // hmm, not used
    )

    ch_versions = ch_versions.mix(SPADES.out.versions)

    MEGAHIT(
        qc_reads.megahit
    )
    
    assembly = SPADES.out.contigs.join(MEGAHIT.out.contigs, remainder: true)
                .collect(sort: metaSorter)
                .map { nothing, contigs, meta -> [meta, contigs] }
    
    ch_versions = ch_versions.mix(MEGAHIT.out.versions)

    // Clean the assembly contigs //
    ASSEMBLY_QC(
        assembly,
        params.reference_genome
    )

    ch_versions = ch_versions.mix(ASSEMBLY_QC.out.versions)

    // Coverage //
    ASSEMBLY_COVERAGE(
        READS_QC.out.qc_reads,
        ASSEMBLY_QC.out.filtered_contigs
    )

    ch_versions = ch_versions.mix(ASSEMBLY_COVERAGE.out.versions)

    // Stats //
    /* The QUAST module was modified to run metaQUAST instead */
    QUAST(
        ASSEMBLY_QC.out.filtered_contigs,
        [ [], [] ], // reference
        [ [], [] ]  // gff
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMiassembler.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMiassembler.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ASSEMBLY_COVERAGE.out.samtools_idxstats.collect{ it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.results.collect { it[1] }.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
