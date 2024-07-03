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
   log.info paramsHelp("nextflow run ebi-metagenomics/miassembler --help")
   exit 0
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.fromPath("$projectDir/assets/mgnify_logo.png")
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
include { FASTQC as FASTQC_BEFORE      } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_AFTER       } from '../modules/nf-core/fastqc/main'
include { MULTIQC                      } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS  } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SPADES                       } from '../modules/nf-core/spades/main'
include { MEGAHIT                      } from '../modules/nf-core/megahit/main'
include { QUAST                        } from '../modules/nf-core/quast/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []


workflow MIASSEMBLER {

    ch_versions = Channel.empty()

    if ( params.samplesheet ) {
        groupReads = { study_accession, reads_accession, fq1, fq2, library_layout, library_strategy, assembly_memory ->
            if (fq2 == []) {
                return tuple(["id": reads_accession,
                              "study_accession": study_accession,
                              "library_strategy": library_strategy,
                              "library_layout": library_layout,
                              "single_end": true,
                              "assembly_memory": assembly_memory ?: params.assembly_memory
                            ],
                            [fq1]
                        )
            } else {
                return tuple(["id": reads_accession,
                              "study_accession": study_accession,
                              "library_strategy": library_strategy,
                              "library_layout": library_layout,
                              "single_end": false,
                              "assembly_memory": assembly_memory ?: params.assembly_memory
                            ],
                            [fq1, fq2])
            }
        }

        samplesheet = Channel.fromList(samplesheetToList(params.samplesheet, "./assets/schema_input.json"))

        // [ study, sample, read1, [read2], library_layout, library_strategy, assembly_memory ]
        fetch_reads_transformed = samplesheet.map(groupReads)

        fetch_reads_transformed.view()

    } else {
        // TODO: remove when the fetch tools get's published on bioconda
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
        fetch_reads_transformed = FETCHTOOL_READS.out.reads.map { meta, reads, library_strategy, library_layout -> {
                [ meta + [
                    //  -- The metadata will be overriden by the parameters -- //
                    "library_strategy": params.library_strategy ?: library_strategy,
                    "library_layout": params.library_layout ?: library_layout,
                    "single_end": params.single_end ?: library_layout == "single"
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

    FASTQC_BEFORE (
        fetch_reads_transformed
    )

    ch_versions = ch_versions.mix(FASTQC_BEFORE.out.versions)

    READS_QC(
        fetch_reads_transformed,
        params.reference_genome
    )

    FASTQC_AFTER (
        READS_QC.out.qc_reads
    )

    /***************************/
    /* Selecting the assembler */
    /***************************/
    /*
        The selection process ensures that:
        - The user selected assembler is always used.
        - Single-end reads are assembled with MEGAHIT, unless specified otherwise.
        - Paired-end reads are assembled with MetaSPAdes, unless specified otherwise
        - An error is raised if the assembler and read layout are incompatible (shouldn't happen...)
    */
    qc_reads_extended = READS_QC.out.qc_reads.map { meta, reads ->
        if ( params.assembler == "megahit" || ( meta.single_end && params.assembler == null ) ) {
            return [ meta + [assembler: "megahit", assembler_version: params.megahit_version], reads]
        } else if ( ["metaspades", "spades"].contains(params.assembler) || ( !meta.single_end && params.assembler == null ) ) {
            def xspades_assembler = params.assembler ?: "metaspades" // Default to "metaspades" if the user didn't select one
            return [ meta + [assembler: xspades_assembler, assembler_version: params.spades_version], reads]
        } else {
            error "Incompatible assembler and/or reads layout. We can't assembly data that is. Reads - single end value: ${meta.single_end}."
        }
    }

    qc_reads_extended.branch { meta, reads ->
        megahit: meta.assembler == "megahit"
        xspades: ["metaspades", "spades"].contains(meta.assembler)
    }.set { qc_reads }

    ch_versions = ch_versions.mix(READS_QC.out.versions)

    /*********************/
    /*     Assembly     */
    /********************/
    /* -- Clarification -- */
    /*
        At the moment, the pipeline only processes one set of reads at a time.
        Therefore, running Spades, metaSpades, or MEGAHIT are mutually exclusive.
        In order to support multiple runs, we need to refactor the code slightly.
        We will need to use `.join()` to keep assemblies together with their corresponding reads.
    */

    SPADES(
        qc_reads.xspades.map { meta, reads -> [meta, reads, [], []] },
        params.assembler ?: "metaspades",
        [], // yml input parameters, which we don't use
        []  // hmm, not used
    )

    ch_versions = ch_versions.mix(SPADES.out.versions)

    MEGAHIT(
        qc_reads.megahit
    )

    assembly = SPADES.out.contigs.mix( MEGAHIT.out.contigs )

    ch_versions = ch_versions.mix(MEGAHIT.out.versions)

    // Clean the assembly contigs //
    ASSEMBLY_QC(
        assembly,
        params.reference_genome
    )

    ch_versions = ch_versions.mix(ASSEMBLY_QC.out.versions)

    // Coverage //
    ASSEMBLY_COVERAGE(
        qc_reads_extended.join( ASSEMBLY_QC.out.filtered_contigs )
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
    if (!params.samplesheet) {
        ch_multiqc_files = ch_multiqc_files.mix(fetch_tool_metadata)
    }
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_BEFORE.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_AFTER.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ASSEMBLY_COVERAGE.out.samtools_idxstats.collect{ it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.results.collect { it[1] }.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        QUAST.out.results.map{meta, _ -> meta}
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
