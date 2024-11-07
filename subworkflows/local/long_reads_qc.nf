include { FASTP as FASTP_LR                      } from '../../modules/nf-core/fastp/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_HUMAN } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_HOST  } from '../../modules/nf-core/minimap2/align/main'

workflow LONG_READS_QC {

    take:
    reads              // [ val(meta), path(reads) ]
    reference_genome   // [ val(meta2), path(reference_genome) ]

    main:
    ch_versions = Channel.empty()

    FASTP_LR(
        reads,
        [],      // no input adapters
        false,   // keep passing reads in the output
        false,   // omit trimmed reads in the output
        false,   // don't merge all reads in the output
        false    // don't trim for polyA
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

    // TODO: add filter if too many reads are removed

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

        MINIMAP2_ALIGN_HUMAN(
            FASTP_LR.out.reads,
            human_reference,
            "human",
            true,    // output bam format
            "bai",   // bam index extension
            false,   // no CIGAR in paf format
            true     // allow for long CIGAR
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

        MINIMAP2_ALIGN_HOST(
            decontaminated_reads,
            host_reference,
            "host",
            true,    // output bam format
            "bai",   // bam index extension
            false,   // no CIGAR in paf format
            true     // allow for long CIGAR
        )

        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_HOST.out.versions)

        decontaminated_reads = MINIMAP2_ALIGN_HOST.out.filtered_fastq
    }

    emit:
    qc_reads = decontaminated_reads
    versions = ch_versions
}
