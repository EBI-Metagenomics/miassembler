include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_HUMAN } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP_ALIGN_HOST  } from '../../modules/nf-core/minimap2/align/main'

workflow LONG_READS_ASSEMBLY_QC {

    take:
    assembly              // [ val(meta), path(assembly_fasta) ]
    host_reference_genome // [ val(meta2), path(host_reference_genome) ] | meta2 contains the name of the reference genome

    main:

    ch_versions = Channel.empty()

    // TODO: post-assembly decontamination could be a separate workflow since it's repeated across different LR tracks
    if ( params.remove_human ) {
        // TODO: double check this after merging SR and LR
        // TODO: check if the id in the mapping below should go back to "human_blast_index_name". If not correct merged version too
        human_reference = Channel.fromPath( "${params.reference_genomes_folder}/${params.human_fasta_prefix}.fna", checkIfExists: true)
            .collect().map {
                files -> [ ["id": params.human_fasta_prefix], files ]
            }

        // TODO: can we change the way human/host are given via prefixes?

        MINIMAP2_ALIGN_HUMAN(
            FLYE.out.fasta,
            human_reference,
            "human",
            true,
            "bai",
            false,
            true
        )

        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_HUMAN.out.versions)

        decontaminated_assembly = MINIMAP2_ALIGN_HUMAN.out.filtered_fastq

    } else {
        decontaminated_assembly = FLYE.out.fasta
    }

    if ( host_reference_genome != null ) {

        host_reference = Channel.fromPath( "${params.reference_genomes_folder}/${host_reference_genome}*", checkIfExists: true)
            .collect().map {
                files -> [ ["id": host_reference_genome], files ]
            }

        MINIMAP_ALIGN_HOST(
            decontaminated_assembly,
            host_reference,
            "host",
            true,
            "bai",
            false,
            true
        )

        ch_versions = ch_versions.mix(MINIMAP_ALIGN_HOST.out.versions)

        decontaminated_assembly = MINIMAP_ALIGN_HOST.out.filtered_fastq
    }

    // temporary just to test the subwf
    emit:
    contigs  = decontaminated_assembly
    versions = ch_versions
}