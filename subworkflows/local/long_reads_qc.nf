include { FASTP_LR                                } from '../../modules/nf-core/fastp/main'
include { RAW_READ_QUALITY_CHECK                  } from '../../modules/local/raw_read_quality_check/'
include { MINIMAP2_ALIGN as HUMAN_DECONTAMINATION } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as HOST_DECONTAMINATION  } from '../../modules/nf-core/minimap2/align/main'

workflow LONG_READS_QC {
    take:
    reads                   // [ val(meta), path(reads) ]
    host_reference_genome   // [ val(meta2), path(reference_genome) ]

    main:
    ch_versions = Channel.empty()

    FASTP_LR(
        reads,
        [],
        false,
        false,
        false,
        false
    )

    ch_versions = ch_versions.mix(FASTP.out.versions)

    RAW_READ_QUALITY_CHECK(
        FASTP.out.json
    )

    decontaminated_reads = channel.empty()

    if ( params.remove_human ) {

        human_reference = Channel.fromPath( "${params.reference_genomes_folder}/${params.human_fasta_prefix}.fna", checkIfExists: true)
            .collect().map {
                files -> [ ["id": params.human_blast_index_name], files ]
            }

        // TODO: can we change the way human/host are given via prefixes?

        HUMAN_DECONTAMINATION(
            FASTP.out.reads,
            human_reference,
            "human",
            true,
            "bai",
            false,
            true
        )

        ch_versions = ch_versions.mix(HUMAN_DECONTAMINATION.out.versions)

        decontaminated_reads = HUMAN_DECONTAMINATION.out.filtered_fastq

    } else {
        decontaminated_reads = FASTP.out.reads
    }

    if ( host_reference_genome != null ) {

        host_reference = Channel.fromPath( "${params.reference_genomes_folder}/${host_reference_genome}*", checkIfExists: true)
            .collect().map {
                files -> [ ["id": host_reference_genome], files ]
            }

        HOST_DECONTAMINATION(
            decontaminated_reads,
            host_reference,
            "host",
            true,
            "bai",
            false,
            true
        )

        ch_versions = ch_versions.mix(HOST_DECONTAMINATION.out.versions)

        decontaminated_reads = HOST_DECONTAMINATION.out.filtered_fastq
    }

    final_reads = decontaminated_reads
                .map{ meta, reads -> {
                        [ meta + [
                            "quality": RAW_READ_QUALITY_CHECK.out.quality.val
                        ], reads ]
                    }
                }

    emit:
    qc_reads = final_reads
    versions = ch_versions
}
