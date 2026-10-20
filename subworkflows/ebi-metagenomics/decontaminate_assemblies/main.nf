include { DECONTAMINATE_CONTIGS as HUMAN_DECONTAMINATE_CONTIGS } from '../decontaminate_contigs/main'
include { DECONTAMINATE_CONTIGS as PHIX_DECONTAMINATE_CONTIGS  } from '../decontaminate_contigs/main'
include { DECONTAMINATE_CONTIGS as HOST_DECONTAMINATE_CONTIGS  } from '../decontaminate_contigs/main'

workflow DECONTAMINATE_ASSEMBLIES {
    /*
    * Microbiome Informatics metagenomics assembly decontamination subworkflow
    *
    * Performs sequential decontamination of assembled contigs against the human genome, PhiX and a genome
    * specified by the user (usually the sample host). It is designed to remove any contigs that come from any
    * of those references.
    *
    * Rationale:
    * 1. PhiX sequences - commonly used as spike-in controls in Illumina sequencing
    * 2. Human sequences - to remove human contamination from samples
    * 3. Host/contaminant sequences - custom contaminants specific to the experimental design
    *
    * Each decontamination step is optional and controlled by metadata parameters, allowing
    * flexible configuration based on sample requirements.
    *
    * Global parameters:
    * - params.reference_genomes_folder: Path to directory containing reference genome files
    *
    * Parameters (via meta map):
    *
    * - meta.phix_reference: Filename prefix of PhiX reference genome (optional)
    *   - If null, PhiX decontamination is skipped
    *   - File must exist in "${params.reference_genomes_folder}/${meta.phix_reference}/${meta.phix_reference}.fna"
    *
    * - meta.human_reference: Filename prefix of human reference genome (optional)
    *   - If null, human decontamination is skipped
    *   - File must exist in "${params.reference_genomes_folder}/${meta.human_reference}/${meta.human_reference}.fna"
    *
    * - meta.contaminant_reference: Filename prefix of custom contaminant reference (optional)
    *   - If null, host/contaminant decontamination is skipped
    *   - File must exist in "${params.reference_genomes_folder}/${meta.contaminant_reference}/${meta.contaminant_reference}.fna"
    *
    * Input:
    * - assembly: Channel [ val(meta), path(assembly_fasta) ]
    *
    * Outputs:
    * - cleaned_contigs: Channel [ val(meta), path(cleaned_assembly_fasta) ]
    * - versions: Channel containing software version information
    *
    */

    take:
    assembly               // [ val(meta), path(assembly_fasta) ]

    main:
    ch_versions = Channel.empty()

    /***************************************************************************/
    /* Perform decontamination from PhiX sequences if requested                */
    /***************************************************************************/
    assembly.branch { meta, _contigs ->
            run_decontamination: meta.phix_reference != null
            skip_decontamination: meta.phix_reference == null
        }
        .set { phix_subdivided_assemblies }

    phix_subdivided_assemblies.run_decontamination
        .map { meta, contigs ->
            [ [meta, contigs], file( "${params.reference_genomes_folder}/${meta.phix_reference}/${meta.phix_reference}.fna", checkIfExists: true )  ]
        }
        .set { ch_phix_decontamination_input }

    PHIX_DECONTAMINATE_CONTIGS(ch_phix_decontamination_input)
    ch_versions = ch_versions.mix(PHIX_DECONTAMINATE_CONTIGS.out.versions)

    phix_cleaned_contigs = phix_subdivided_assemblies.skip_decontamination.mix(
        PHIX_DECONTAMINATE_CONTIGS.out.cleaned_contigs
    )

    /***************************************************************************/
    /* Perform decontamination from human sequences if requested               */
    /***************************************************************************/
    phix_cleaned_contigs
        .branch { meta, _contigs ->
            run_decontamination: meta.human_reference != null
            skip_decontamination: meta.human_reference == null
        }
        .set { human_subdivided_assemblies }

    human_subdivided_assemblies.run_decontamination
        .map { meta, contigs ->
            [ [meta, contigs], file( "${params.reference_genomes_folder}/${meta.human_reference}/${meta.human_reference}.fna", checkIfExists: true ) ]
        }
        .set { ch_human_decontamination_input }

    HUMAN_DECONTAMINATE_CONTIGS(ch_human_decontamination_input)
    ch_versions = ch_versions.mix(HUMAN_DECONTAMINATE_CONTIGS.out.versions)

    human_cleaned_contigs = human_subdivided_assemblies.skip_decontamination.mix(
        HUMAN_DECONTAMINATE_CONTIGS.out.cleaned_contigs
    )

    /***************************************************************************/
    /* Perform decontamination from arbitrary contaminant sequences            */
    /***************************************************************************/
    human_cleaned_contigs
        .branch { meta, _contigs ->
            run_decontamination: meta.contaminant_reference != null
            skip_decontamination: meta.contaminant_reference == null
        }
        .set { subdivided_assemblies }

    subdivided_assemblies.run_decontamination
        .map { meta, contigs ->
            [ [meta, contigs], file( "${params.reference_genomes_folder}/${meta.contaminant_reference}/${meta.contaminant_reference}.fna", checkIfExists: true ) ]
        }
        .set { ch_decontamination_input }

    HOST_DECONTAMINATE_CONTIGS(ch_decontamination_input)
    ch_versions = ch_versions.mix(HOST_DECONTAMINATE_CONTIGS.out.versions)

    cleaned_contigs = subdivided_assemblies.skip_decontamination.mix(
        HOST_DECONTAMINATE_CONTIGS.out.cleaned_contigs
    )

    emit:
    cleaned_contigs                    = cleaned_contigs
    phix_contaminated_contigs_tsv      = PHIX_DECONTAMINATE_CONTIGS.out.mapped_contigs_tsv
    human_contaminated_contigs_tsv     = HUMAN_DECONTAMINATE_CONTIGS.out.mapped_contigs_tsv
    host_contaminated_contigs_tsv      = HOST_DECONTAMINATE_CONTIGS.out.mapped_contigs_tsv
    phix_contaminated_contigs_tsv_mqc  = PHIX_DECONTAMINATE_CONTIGS.out.mapped_contigs_tsv_mqc
    human_contaminated_contigs_tsv_mqc = HUMAN_DECONTAMINATE_CONTIGS.out.mapped_contigs_tsv_mqc
    host_contaminated_contigs_tsv_mqc  = HOST_DECONTAMINATE_CONTIGS.out.mapped_contigs_tsv_mqc
    versions                           = ch_versions
}
