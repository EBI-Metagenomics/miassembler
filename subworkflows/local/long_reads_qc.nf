include { FASTP as FASTP_LR                      } from '../../modules/nf-core/fastp/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_HUMAN } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_HOST  } from '../../modules/nf-core/minimap2/align/main'

workflow LONG_READS_QC {

    take:
    input_reads        // [ val(meta), path(reads) ]

    main:
    ch_versions = Channel.empty()

    FASTP_LR(
        input_reads,
        [],      // no input adapters
        false,   // keep passing reads in the output
        false,   // omit trimmed reads in the output
        false,   // don't merge all reads in the output
        false    // don't trim for polyA
    )

    ch_versions = ch_versions.mix(FASTP_LR.out.versions)

    reads_json = FASTP_LR.out.reads.join( FASTP_LR.out.json )

    reads_quality_levels = reads_json.map { meta, reads, json ->
        def json_txt = new groovy.json.JsonSlurper().parseText(json.text)

        def q20_percentage = json_txt?.summary?.before_filtering?.q20_rate ?: 0;

        if ( q20_percentage >= params.long_reads_pacbio_quality_threshold ) {
            return [ meta + [quality: "high"], reads]
        } else {
            return [ meta + [quality: "low"], reads]
        }
    }

    // TODO: add filter if too many reads are removed

    /***************************************************************************/
    /* Perform decontamination from human sequences if requested               */
    /***************************************************************************/
    reads_quality_levels
        .branch { meta, _reads ->
            run_decontamination: meta.human_reference != null
            skip_decontamination: meta.human_reference == null
        }
        .set { human_subdivided_reads }

    human_subdivided_reads.run_decontamination
        .multiMap { meta, reads ->
            reads: [meta, reads]
            reference: [ [id:"human"], "${params.reference_genomes_folder}/${meta.human_reference}" ]
        }
        .set { ch_human_decontamination_input }

    MINIMAP2_ALIGN_HUMAN(
        ch_human_decontamination_input.reads,
        ch_human_decontamination_input.reference,
        "human",
        "fastq", // out sequence extension
        true,    // output bam format
        "bai",   // bam index extension
        false,   // no CIGAR in paf format
        true     // allow for long CIGAR
    )

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_HUMAN.out.versions)

    human_cleaned_reads = human_subdivided_reads.skip_decontamination.mix(
        MINIMAP2_ALIGN_HUMAN.out.filtered_output
    )

    /***************************************************************************/
    /* Perform decontamination from arbitrary contaminant sequences            */
    /***************************************************************************/
    human_cleaned_reads
        .branch { meta, _reads ->
            run_decontamination: meta.contaminant_reference != null
            skip_decontamination: meta.contaminant_reference == null
        }
        .set { subdivided_reads }

    subdivided_reads.run_decontamination
        .multiMap { meta, reads ->
            reads: [meta, reads]
            reference: [ [id:file(meta.contaminant_reference).baseName], "${params.reference_genomes_folder}/${meta.contaminant_reference}" ]
        }
        .set { ch_decontamination_input }

    MINIMAP2_ALIGN_HOST(
        ch_decontamination_input.reads,
        ch_decontamination_input.reference,
        "host",
        "fastq", // out sequence extension
        true,    // output bam format
        "bai",   // bam index extension
        false,   // no CIGAR in paf format
        true     // allow for long CIGAR
    )

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_HOST.out.versions)

    decontaminated_reads = subdivided_reads.skip_decontamination.mix(
        MINIMAP2_ALIGN_HOST.out.filtered_output
    )

    emit:
    qc_reads   = decontaminated_reads
    fastp_json = FASTP_LR.out.json
    versions   = ch_versions
}
