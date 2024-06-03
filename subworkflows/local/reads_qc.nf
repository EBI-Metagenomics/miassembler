include { FASTP                                             } from '../../modules/nf-core/fastp/main'
include { BWAMEM2DECONTNOBAMS as HUMAN_PHIX_DECONTAMINATION } from '../../modules/ebi-metagenomics/bwamem2decontnobams/main'
include { BWAMEM2DECONTNOBAMS as HOST_DECONTAMINATION       } from '../../modules/ebi-metagenomics/bwamem2decontnobams/main'

workflow READS_QC {

    take:
    reads                 // [ val(meta), path(reads) ]
    host_reference_genome // [ val(meta2), path(reference_genome) ] | meta2 contains the name of the reference genome

    main:
    ch_versions = Channel.empty()

    FASTP(
        reads,
        [],
        false,
        false,
        false
    )

    ch_versions = ch_versions.mix(FASTP.out.versions)

    decontaminated_reads = channel.empty()

    if ( params.remove_human_phix ) {

        ch_bwamem2_human_phix_refs = Channel.fromPath( "${params.bwamem2_reference_genomes_folder}/${params.human_phix_blast_index_name}*", checkIfExists: true)
            .collect().map {
                files -> [ ["id": params.human_phix_blast_index_name], files ]
            }

        HUMAN_PHIX_DECONTAMINATION(
            FASTP.out.reads,
            ch_bwamem2_human_phix_refs
        )

        ch_versions = ch_versions.mix(HUMAN_PHIX_DECONTAMINATION.out.versions)

        decontaminated_reads = HUMAN_PHIX_DECONTAMINATION.out.decont_reads

    } else {
        decontaminated_reads = FASTP.out.reads
    }

    if ( host_reference_genome != null ) {

        ch_bwamem2_host_refs = Channel.fromPath( "${params.bwamem2_reference_genomes_folder}/${host_reference_genome}*", checkIfExists: true)
            .collect().map {
                files -> [ ["id": host_reference_genome], files ]
            }

        HOST_DECONTAMINATION(
            HUMAN_PHIX_DECONTAMINATION.out.decont_reads,
            ch_bwamem2_host_refs
        )

        ch_versions = ch_versions.mix(HOST_DECONTAMINATION.out.versions)

        decontaminated_reads = HOST_DECONTAMINATION.out.decont_reads
    }

    emit:
    qc_reads = decontaminated_reads
    versions = ch_versions
}
