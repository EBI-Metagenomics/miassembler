include { DECONTAMINATE_CONTIGS as HUMAN_DECONTAMINATE_CONTIGS } from '../../subworkflows/ebi-metagenomics/decontaminate_contigs/main'
include { DECONTAMINATE_CONTIGS as HOST_DECONTAMINATE_CONTIGS  } from '../../subworkflows/ebi-metagenomics/decontaminate_contigs/main'

workflow LONG_READS_ASSEMBLY_QC {

    take:
    assembly            // [ val(meta), path(assembly_fasta) ]
    reference_genome    // [ val(reference_file) ]

    main:

    ch_versions = Channel.empty()

    if ( params.remove_human ) {

        ch_human_ref = Channel
            .fromPath(
                "${params.reference_genomes_folder}/${params.human_fasta_prefix}.f*a",
                checkIfExists: true
            )
            .collect() // TODO I don't think collect() is needed here, but it was used in the original code
            .map {
                file -> [ ["id": params.human_fasta_prefix], file ]
            }

        HUMAN_DECONTAMINATE_CONTIGS(
            assembly,
            ch_human_ref
        )

        ch_versions = ch_versions.mix(HUMAN_DECONTAMINATE_CONTIGS.out.versions)

        decontaminated_assembly = HUMAN_DECONTAMINATE_CONTIGS.out.cleaned_contigs

    } else {
        decontaminated_assembly = assembly
    }

    if ( reference_genome ) {

        ch_contaminat_ref = Channel
            .fromPath(
                "${params.reference_genomes_folder}/${reference_genome}",
                checkIfExists: true
            )
            .map {
                file -> [ ["id": reference_genome], file ]
            }

        HOST_DECONTAMINATE_CONTIGS(
            decontaminated_assembly,
            ch_contaminat_ref
        )

        ch_versions = ch_versions.mix(HOST_DECONTAMINATE_CONTIGS.out.versions)

        decontaminated_assembly = HOST_DECONTAMINATE_CONTIGS.out.cleaned_contigs
    }

    emit:
    contigs  = decontaminated_assembly
    versions = ch_versions
}