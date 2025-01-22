include { FASTP as FASTP_LR                      } from '../../modules/nf-core/fastp/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_HUMAN } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_HOST  } from '../../modules/nf-core/minimap2/align/main'

workflow LONG_READS_QC {

    take:
    input_reads        // [ val(meta), path(reads) ]
    reference_genome   // [ val(meta2), path(reference_genome) ]

    main:
    def ch_versions = Channel.empty()
    def human_reference = Channel.empty()
    def host_reference = Channel.empty()

    FASTP_LR(
        input_reads,
        [],      // no input adapters
        false,   // keep passing reads in the output
        false,   // omit trimmed reads in the output
        false,   // don't merge all reads in the output
        false    // don't trim for polyA
    )

    ch_versions = ch_versions.mix(FASTP_LR.out.versions)

    def reads_json = FASTP_LR.out.reads.join( FASTP_LR.out.json )

    def reads_quality_levels = reads_json.map { meta, reads, json ->
        def json_txt = new groovy.json.JsonSlurper().parseText(json.text)
        
        def q20_percentage = json_txt?.summary?.before_filtering?.q20_rate ?: 0;

        if ( q20_percentage >= params.long_reads_pacbio_quality_threshold ) {
            return [ meta + [quality: "high"], reads]
        } else {
            return [ meta + [quality: "low"], reads]
        }
    }

    // TODO: add filter if too many reads are removed

    def decontaminated_reads = channel.empty()

    if ( params.remove_human ) {
        // TODO: make this consistent with short_reads
        // can we use the same flag, even if one has phix but not the other?
        // Check file extensions too

        human_reference = Channel.fromPath(
            "${params.reference_genomes_folder}/${params.human_fasta_prefix}.f*a", checkIfExists: true)
            .collect().map {
                files -> [ ["id": params.human_fasta_prefix], files ]
            }

        // TODO: can we change the way human/host are given via prefixes?

        MINIMAP2_ALIGN_HUMAN(
            reads_quality_levels,
            human_reference,
            "human",
            "fastq", // out sequence extension
            true,    // output bam format
            "bai",   // bam index extension
            false,   // no CIGAR in paf format
            true     // allow for long CIGAR
        )

        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_HUMAN.out.versions)

        decontaminated_reads = MINIMAP2_ALIGN_HUMAN.out.filtered_output

    } else {
        decontaminated_reads = reads_quality_levels
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
            "fastq", // out sequence extension
            true,    // output bam format
            "bai",   // bam index extension
            false,   // no CIGAR in paf format
            true     // allow for long CIGAR
        )

        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_HOST.out.versions)

        decontaminated_reads = MINIMAP2_ALIGN_HOST.out.filtered_output
    }

    emit:
    qc_reads = decontaminated_reads
    versions = ch_versions
}
