/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT PLUGINS AND OTHER BITS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsSummaryLog; paramsSummaryMap; samplesheetToList; paramsHelp } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { MULTIQC as MULTIQC_STUDY    } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_RUN      } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT THE MAIN ENTRY POINT WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOWS
//
include { SHORT_READS_ASSEMBLER       } from '../workflows/short_reads_assembler'
include { LONG_READS_ASSEMBLER        } from '../workflows/long_reads_assembler'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DOWNLOAD_FROM_FIRE as DOWNLOAD_FROM_FIRE_SHORT_READS } from '../modules/ebi-metagenomics/downloadfromfire/main'
include { DOWNLOAD_FROM_FIRE as DOWNLOAD_FROM_FIRE_LONG_READS  } from '../modules/ebi-metagenomics/downloadfromfire/main'
include { FETCHTOOL_READS                                      } from '../modules/local/fetchtool_reads'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MIASSEMBLER {

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PRINT PARAMS SUMMARY
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
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

    def ch_multiqc_config          = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    def ch_multiqc_custom_config   = params.multiqc_config ? file( params.multiqc_config, checkIfExists: true ) : []
    def ch_multiqc_logo            = params.multiqc_logo   ? file( params.multiqc_logo, checkIfExists: true ) : file("$projectDir/assets/mgnify_logo.png", checkIfExists: true)
    def ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    def ch_versions = Channel.empty()
    def fetch_tool_metadata = Channel.empty()
    def fetch_reads_transformed = Channel.empty()

    // Print parameter summary log to screen
    log.info(logo + paramsSummaryLog(workflow) + citation)

    // ***************************************************************************** //
    // Custom validation for human decontamination of reads and assembled contigs   //
    // ***************************************************************************** //

    if (params.skip_human_decontamination) {
        log.warn("Human sequences will not be removed from raw reads or assembled contigs.")

        if (params.human_reference) {
            error "Setting 'skip_human_decontamination = true' is incompatible with providing 'params.human_reference'."
        }
    }

    if (params.samplesheet) {
        def groupReads = { study_accession, reads_accession, fq1, fq2, library_layout, library_strategy, platform, assembler, assembly_memory, assembler_config, contaminant_reference, human_reference, phix_reference, lambdaphage_reference ->

            def human_reference_path = human_reference ?: params.human_reference
            if (!params.skip_human_decontamination && human_reference_path == null) {
                error "Invalid row, skip_human_decontamination is false but there is no human reference on row: ${study_accession}, ${reads_accession}."
            }

            if (fq2 == []) {
                return tuple(
                    [
                        "id": reads_accession,
                        "study_accession": study_accession,
                        "single_end": true,
                        "library_layout": library_layout,
                        "library_strategy": library_strategy,
                        "platform": params.platform ?: platform,
                        "assembler": assembler ?: params.assembler,
                        "assembly_memory": assembly_memory ?: params.assembly_memory,
                        "assembler_config": assembler_config ?: params.long_reads_assembler_config,
                        "contaminant_reference": contaminant_reference ?: params.contaminant_reference,
                        "human_reference": human_reference_path, // -> if this value is null (which is not the same as an empty string) the decontamination won't be executed
                        "phix_reference": phix_reference ?: params.phix_reference,
                        "lambdaphage_reference": lambdaphage_reference ?: params.lambdaphage_reference
                    ],
                    [fq1]
                )
            } else {
                return tuple(
                    [
                        "id": reads_accession,
                        "study_accession": study_accession,
                        "single_end": false,
                        "library_layout": library_layout,
                        "library_strategy": library_strategy,
                        "platform": params.platform ?: platform,
                        "assembler": assembler ?: params.assembler,
                        "assembly_memory": assembly_memory ?: params.assembly_memory,
                        "assembler_config": assembler_config ?: params.long_reads_assembler_config,
                        "contaminant_reference": contaminant_reference ?: params.contaminant_reference,
                        "human_reference": human_reference_path, // -> if this value is null (which is not the same as an empty string) the decontamination won't be executed
                        "phix_reference": phix_reference ?: params.phix_reference,
                        "lambdaphage_reference": lambdaphage_reference ?: params.lambdaphage_reference
                    ],
                    [fq1, fq2]
                )
            }
        }

        def samplesheet = Channel.fromList(samplesheetToList(params.samplesheet, "./assets/schema_input.json"))

        // [ study, sample, read1, [read2], library_layout, library_strategy, platform, assembly_memory]
        fetch_reads_transformed = samplesheet.map(groupReads)
    }
    else {
        // TODO: remove when the fetch tools get's published on bioconda
        def fetch_tool_config = file("${projectDir}/assets/fetch_tool_anonymous.json", checkIfExists: true)

        if (params.private_study) {
            fetch_tool_config = file("${projectDir}/assets/fetch_tool_credentials.json", checkIfExists: true)
        }

        FETCHTOOL_READS(
            [[id: params.reads_accession], params.study_accession, params.reads_accession],
            fetch_tool_config
        )

        ch_versions = ch_versions.mix(FETCHTOOL_READS.out.versions)

        if (!params.skip_human_decontamination && params.human_reference == null) {
                error "Human decontamination is enabled but no human reference is provided. Please specify 'human_reference' parameter or set 'skip_human_decontamination = true'."
            }

        // Push the library strategy into the meta of the reads, this is to make it easier to handle downstream
        fetch_reads_transformed = FETCHTOOL_READS.out.reads.map { meta, reads, library_strategy, library_layout, platform ->
            {
                [
                    meta + [
                        "assembler": params.assembler,
                        "assembler_config": params.long_reads_assembler_config,
                        "assembly_memory": params.assembly_memory,
                        "library_strategy": params.library_strategy ?: library_strategy,
                        "library_layout": params.library_layout ?: library_layout,
                        "single_end": params.single_end ?: library_layout == "single",
                        "platform": params.platform ?: platform,
                        "contaminant_reference": params.contaminant_reference,
                        "human_reference": params.skip_human_decontamination ? null : params.human_reference,
                        "phix_reference": params.phix_reference,
                        "lambdaphage_reference": params.lambdaphage_reference
                    ],
                    reads
                ]
            }
        }

        // Metadata for MultiQC
        fetch_tool_metadata = FETCHTOOL_READS.out.metadata_tsv
            .map { it[1] }
            .collectFile(
                name: 'fetch_tool_mqc.tsv',
                newLine: true,
                keepHeader: true,
                skip: 1
            )
    }

    /*******************************************/
    /* Selecting the assembly pipeline flavour */
    /*******************************************/

    def classified_reads = fetch_reads_transformed.map { meta, reads ->
        // Long reads //
        if ( ["ont", "pb"].contains( meta.platform ) ) {
            return [ meta + [long_reads: true], reads]
        // Short reads //
        } else {
            return [ meta + [short_reads: true], reads]
        }
    }

    classified_reads
        .branch { meta, _reads ->
            short_reads: meta.short_reads
            long_reads: meta.long_reads
        }
    .set { reads_to_assemble }

    /**********************/
    // DOWNLOAD THE READS //
    /**********************/

    def short_reads = reads_to_assemble.short_reads
    def long_reads = reads_to_assemble.long_reads

    // If running on EBI infrastructure, and on samplehseets otherwise the fetch tool will kick in //
    if (params.samplesheet && params.use_fire_download) {
        /*
         * For private studies we need to bypass Nextflow S3 integration until https://github.com/nextflow-io/nextflow/issues/4873 is fixed
         * The EBI parameter is needed as this only works on EBI network, FIRE is not accessible otherwise
        */
        DOWNLOAD_FROM_FIRE_SHORT_READS(
            short_reads
        )
        ch_versions = ch_versions.mix(DOWNLOAD_FROM_FIRE_SHORT_READS.out.versions.first())

        short_reads = DOWNLOAD_FROM_FIRE_SHORT_READS.out.downloaded_files

        DOWNLOAD_FROM_FIRE_LONG_READS(
            long_reads
        )
        ch_versions = ch_versions.mix(DOWNLOAD_FROM_FIRE_LONG_READS.out.versions.first())

        long_reads = DOWNLOAD_FROM_FIRE_LONG_READS.out.downloaded_files
    }

    /***************************************/
    /* Assemble short reads and long reads */
    /***************************************/

    SHORT_READS_ASSEMBLER(
        short_reads
    )

    ch_versions = ch_versions.mix(SHORT_READS_ASSEMBLER.out.versions)

    LONG_READS_ASSEMBLER(
        long_reads
    )

    ch_versions = ch_versions.mix(LONG_READS_ASSEMBLER.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    def workflow_summary    = WorkflowMiassembler.paramsSummaryMultiqc(workflow, summary_params)
    def ch_workflow_summary = Channel.value(workflow_summary)

    def methods_description    = WorkflowMiassembler.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    def ch_methods_description = Channel.value(methods_description)

    def ch_multiqc_base_files = Channel.empty()
    ch_multiqc_base_files = ch_multiqc_base_files.mix( CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect() )
    ch_multiqc_base_files = ch_multiqc_base_files.mix( ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml') )
    ch_multiqc_base_files = ch_multiqc_base_files.mix( ch_methods_description.collectFile(name: 'methods_description_mqc.yaml') )

    if (params.skip_human_decontamination) {
        ch_multiqc_base_files = ch_multiqc_base_files.mix( channel.from( file("$projectDir/assets/human_decontamination_mqc.html", checkIfExists: true) ) )
    }

    /**************************************/
    /* MultiQC report for the whole study */
    /**************************************/

    def meta_by_study = { meta, result_artifact ->
        [meta.subMap("study_accession"), result_artifact]
    }

    // Helper method for the MultiQC aggregation by study and runs //
    def combineFiles = { meta, fastqc_before, fastqc_after, assembly_coverage, quast ->
        // Flatten the fastqc_before and fastqc_after lists
        def flattened_fastqc_before = fastqc_before instanceof List ? fastqc_before.flatten() : [fastqc_before]
        def flattened_fastqc_after = fastqc_after instanceof List ? fastqc_after.flatten() : [fastqc_after]

        // Combine all elements into a single list
        def all_files = flattened_fastqc_before + flattened_fastqc_after
        if (assembly_coverage) {
            all_files += [assembly_coverage]
        }
        if (quast) {
            all_files += [quast]
        }
        // Produce a tuple with meta and the flattened list of files
        return all_files.flatten().collect { file ->
            [meta, file]
        }
    }

    def ch_multiqc_study_tools_files = Channel.empty()

    def fastqc_before_zip = SHORT_READS_ASSEMBLER.out.fastqc_before_zip
        .mix(LONG_READS_ASSEMBLER.out.fastqc_before_zip)
    def fastqc_after_zip = SHORT_READS_ASSEMBLER.out.fastqc_after_zip
        .mix(LONG_READS_ASSEMBLER.out.fastqc_after_zip)
    def assembly_coverage_samtools_idxstats = SHORT_READS_ASSEMBLER.out.assembly_coverage_samtools_idxstats
        .mix(LONG_READS_ASSEMBLER.out.assembly_coverage_samtools_idxstats)
    def quast_results = SHORT_READS_ASSEMBLER.out.quast_results
        .mix(LONG_READS_ASSEMBLER.out.quast_results)

    def study_multiqc_files = fastqc_before_zip.map(meta_by_study)
        .join(fastqc_after_zip.map(meta_by_study))
        .join(assembly_coverage_samtools_idxstats.map(meta_by_study), remainder: true) // the assembly step could fail
        .join(quast_results.map(meta_by_study), remainder: true)                       // the assembly step could fail

    ch_multiqc_study_tools_files = study_multiqc_files.flatMap(combineFiles).groupTuple()

    // TODO: add the fetch tool log file

    MULTIQC_STUDY(
        ch_multiqc_base_files.collect(),
        ch_multiqc_study_tools_files,
        ch_multiqc_config,
        ch_multiqc_custom_config,
        ch_multiqc_logo,
        [],
        []
    )

    /**************************/
    /* MultiQC report per run */
    /*************************/

    def meta_by_run = { meta, result_artifact ->
        [meta.subMap("study_accession", "id", "assembler", "assembler_version"), result_artifact]
    }

    def run_multiqc_files = SHORT_READS_ASSEMBLER.out.fastqc_before_zip.map(meta_by_run)
        .join(SHORT_READS_ASSEMBLER.out.fastqc_after_zip.map(meta_by_run))
        .join(SHORT_READS_ASSEMBLER.out.assembly_coverage_samtools_idxstats.map(meta_by_run), remainder: true) // the assembly step could fail
        .join(SHORT_READS_ASSEMBLER.out.quast_results.map(meta_by_run), remainder: true)                       // the assembly step could fail

    // Filter out the non-assembled runs //
    def ch_multiqc_run_tools_files = run_multiqc_files.filter { _meta, _fastqc_before, _fastqc_after, assembly_coverage, quast -> {
            return assembly_coverage != null && quast != null
        }
    }.flatMap(combineFiles).groupTuple()

    // TODO: add the fetch tool log file
    MULTIQC_RUN(
        ch_multiqc_base_files.collect(),
        ch_multiqc_run_tools_files,
        ch_multiqc_config,
        ch_multiqc_custom_config,
        ch_multiqc_logo,
        [],
        []
    )

    /*****************************/
    /* End of execution reports */
    /****************************/

    // TODO: we need to add LR end-of-run reports

    // Short reads asssembled runs //
    SHORT_READS_ASSEMBLER.out.assembly_coverage_samtools_idxstats
        .map { meta, __ ->
            {
                return "${meta.id},${meta.assembler},${meta.assembler_version}"
            }
        }
        .collectFile(name: "assembled_runs.csv", storeDir: "${params.outdir}", newLine: true, cache: false)

    // Short reads and assembly QC failed //

    def short_reads_qc_failed_entries = SHORT_READS_ASSEMBLER.out.qc_failed_all.map {
        meta, __ -> {
            if (meta.low_reads_count) {
                return "${meta.id},low_reads_count"
            }
            if (meta.filter_ratio_threshold_exceeded) {
                return "${meta.id},filter_ratio_threshold_exceeded"
            }
            if (meta.too_few_contigs) {
                return "${meta.id},too_few_contigs"
            }
            error("Unexpected. meta: ${meta}")
        }
    }

    short_reads_qc_failed_entries.collectFile(name: "qc_failed_runs.csv", storeDir: "${params.outdir}", newLine: true, cache: false)
}
