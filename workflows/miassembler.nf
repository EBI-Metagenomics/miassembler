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
include { FETCHTOOL_READS   } from '../modules/local/fetchtool_reads'
include { PRE_ASSEMBLY_QC   } from '../subworkflows/local/pre_assembly_qc'
include { CLEAN_ASSEMBLY    } from '../subworkflows/local/clean_assembly'
include { ASSEMBLY_COVERAGE } from '../subworkflows/local/assembly_coverage'

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

workflow MIASSEMBLER {

    ch_versions = Channel.empty()

    // Download reads //
    FETCHTOOL_READS(
        [ [id: params.reads_accession], params.study_accession, params.reads_accession ],
        file("$projectDir/assets/fetch_tool_anonymous.json")
    )

    ch_versions = ch_versions.mix(FETCHTOOL_READS.out.versions)

    FASTQC (
        FETCHTOOL_READS.out.reads
    )

    ch_versions = ch_versions.mix(FASTQC.out.versions)
    
    // standard human+phiX genomes dbs
    ch_bwa_humanPhiX_refs = Channel.fromPath( "$params.bwa_reference_genomes_folder/" + params.default_reference_genome + "*", 
        checkIfExists: true).collect().map {
            files -> [ ["id": params.default_reference_genome], files ]
        }
    ch_blast_humanPhiX_refs = Channel.fromPath( "$params.blast_reference_genomes_folder/" + params.default_reference_genome + "*", 
        checkIfExists: true).collect().map {
            files -> [ ["id": params.default_reference_genome], files ]
        }
    
    // Perform QC on reads //
    PRE_ASSEMBLY_QC(
        FETCHTOOL_READS.out.reads, 
        ch_bwa_humanPhiX_refs,
        // ch_bwa_host_refs
        params.reference_genome
    )
    
    ch_versions = ch_versions.mix(PRE_ASSEMBLY_QC.out.versions)

    // Assembly //
    assembly = Channel.empty()

    if ( params.assembler == "metaspades" || params.assembler == "spades" ) {

        SPADES(
            PRE_ASSEMBLY_QC.out.cleaned_reads.map { meta, reads -> [meta, reads, [], []] },
            params.assembler,
            [], // yml input parameters, which we don't use
            []  // hmm, not used
        )

        assembly = SPADES.out.contigs
        ch_versions = ch_versions.mix(SPADES.out.versions)

    } else if ( params.assembler == "megahit" ) {

        MEGAHIT(
            PRE_ASSEMBLY_QC.out.cleaned_reads
        )

        assembly = MEGAHIT.out.contigs
        ch_versions = ch_versions.mix(MEGAHIT.out.versions)

    } else {
        // TODO: raise ERROR, it shouldn't happen as the options are validated by nf-validation
    }

    // Clean the assembly contigs //
    CLEAN_ASSEMBLY(
        assembly,
        ch_blast_humanPhiX_refs,
        // ch_blast_host_refs
        params.reference_genome
    )

    ch_versions = ch_versions.mix(CLEAN_ASSEMBLY.out.versions)

    // Coverage //
    ASSEMBLY_COVERAGE(
        PRE_ASSEMBLY_QC.out.cleaned_reads,
        CLEAN_ASSEMBLY.out.filtered_contigs
    )

    ch_versions = ch_versions.mix(ASSEMBLY_COVERAGE.out.versions)

    // Stats //
    /* The QUAST module was modified to run metaQUAST instead */
    QUAST(
        CLEAN_ASSEMBLY.out.filtered_contigs,
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
