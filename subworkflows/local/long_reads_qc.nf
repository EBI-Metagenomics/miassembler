include { FASTP as FASTP_LR                       } from '../../modules/nf-core/fastp/main'
include { MINIMAP2_ALIGN as HUMAN_DECONTAMINATION } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as HOST_DECONTAMINATION  } from '../../modules/nf-core/minimap2/align/main'

workflow LONG_READS_QC {

    take:
    reads                   // [ val(meta), path(reads) ]
    reference_genome   // [ val(meta2), path(reference_genome) ]

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

    ch_versions = ch_versions.mix(FASTP_LR.out.versions)

    quality_levels_ch = FASTP_LR.out.json.map { meta, json -> {
        json_txt = new JsonSlurper().parseText(json.text)
        q20bases = json_txt?.summary?.before_filtering?.q20_bases ?: 0;
        total_bases = json_txt?.summary?.before_filtering?.total_bases ?: 0;

        q20_percentage = q20_bases / total_bases * 100

        quality = [
            "high_quality": q20_percentage >= 80, 
            "low_quality": q20_percentage < 80,
        ]
        return [meta, quality]
        } 
    }

    RAW_READ_QUALITY_CHECK(
        FASTP_LR.out.json
    )

    decontaminated_reads = channel.empty()

    if ( params.remove_human ) {
        // TODO: make this consistent with short_reads
        // can we use the same flag, even if one has phix but not the other?
        // Check file extensions too

        human_reference = Channel.fromPath( "${params.reference_genomes_folder}/${params.human_fasta_prefix}.fna", checkIfExists: true)
            .collect().map {
                files -> [ ["id": params.human_fasta_prefix], files ]
            }

        // TODO: can we change the way human/host are given via prefixes?

        HUMAN_DECONTAMINATION(
            FASTP_LR.out.reads,
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
        decontaminated_reads = FASTP_LR.out.reads
    }

    if ( reference_genome != null ) {

        host_reference = Channel.fromPath( "${params.reference_genomes_folder}/${reference_genome}*", checkIfExists: true)
            .collect().map {
                files -> [ ["id": reference_genome], files ]
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
