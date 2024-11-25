include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_HUMAN } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP_ALIGN_HOST   } from '../../modules/nf-core/minimap2/align/main'

workflow LONG_READS_ASSEMBLY_QC {

    take:
    assembly            // [ val(meta), path(assembly_fasta) ]
    reference_genome    // [ val(meta2), path(reference_genome) ] | meta2 contains the name of the reference genome

    main:

    ch_versions = Channel.empty()
    decontaminated_assembly = Channel.empty()

    if ( params.remove_human ) {
        human_reference = Channel.fromPath(
            "${params.reference_genomes_folder}/${params.human_fasta_prefix}.fna", checkIfExists: true)
            .collect().map {
                files -> [ ["id": params.human_fasta_prefix], files ]
            }

        MINIMAP2_ALIGN_HUMAN(
            assembly,
            human_reference,
            "human_post",
            "fasta", // alignment output extension
            false    // no CIGAR in paf format
        )

        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_HUMAN.out.versions)

        decontaminated_assembly = MINIMAP2_ALIGN_HUMAN.out.filtered_output

    } else {
        decontaminated_assembly = assembly
    }

    if ( reference_genome != null ) {

        host_reference = Channel.fromPath( "${params.reference_genomes_folder}/${reference_genome}*", checkIfExists: true)
            .collect().map {
                files -> [ ["id": reference_genome], files ]
            }

        MINIMAP_ALIGN_HOST(
            decontaminated_assembly,
            human_reference,
            "host_post",
            "fasta", // alignment output extension
            false    // no CIGAR in paf format
        )

        ch_versions = ch_versions.mix(MINIMAP_ALIGN_HOST.out.versions)

        decontaminated_assembly = MINIMAP_ALIGN_HOST.out.filtered_output
    }

    emit:
    contigs  = decontaminated_assembly
    versions = ch_versions
}