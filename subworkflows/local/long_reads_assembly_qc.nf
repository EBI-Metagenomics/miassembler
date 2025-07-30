include { DECONTAMINATE_CONTIGS as HUMAN_DECONTAMINATE_CONTIGS } from '../../subworkflows/ebi-metagenomics/decontaminate_contigs/main'
include { DECONTAMINATE_CONTIGS as HOST_DECONTAMINATE_CONTIGS  } from '../../subworkflows/ebi-metagenomics/decontaminate_contigs/main'

workflow LONG_READS_ASSEMBLY_QC {

    take:
    assembly            // [ val(meta), path(assembly_fasta) ]

    main:
    ch_versions = Channel.empty()

    /***************************************************************************/
    /* Perform decontamination from human sequences if requested               */
    /***************************************************************************/
    assembly
        .branch { meta, _contigs ->
            run_decontamination: meta.human_reference != null
            skip_decontamination: meta.human_reference == null
        }
        .set { human_subdivided_assemblies }

    human_subdivided_assemblies.run_decontamination
        .map { meta, contigs ->
            [ [meta, contigs], file("${params.reference_genomes_folder}/${meta.human_reference}/minimap2/${meta.human_reference}.fna.mmi", checkIfExists: true) ]
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
            [ [meta, contigs], file("${params.reference_genomes_folder}/${meta.contaminant_reference}/minimap2/${meta.contaminant_reference}.fna.mmi", checkIfExists: true) ]
        }
        .set { ch_decontamination_input }

    HOST_DECONTAMINATE_CONTIGS(ch_decontamination_input)

    cleaned_contigs = subdivided_assemblies.skip_decontamination.mix(
        HOST_DECONTAMINATE_CONTIGS.out.cleaned_contigs
    )

    ch_versions = ch_versions.mix(HOST_DECONTAMINATE_CONTIGS.out.versions)

    emit:
    contigs  = cleaned_contigs
    versions = ch_versions
}
