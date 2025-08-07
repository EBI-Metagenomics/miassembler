include { FASTP as FASTP_LR                      } from '../../modules/nf-core/fastp/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_HUMAN } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_HOST  } from '../../modules/nf-core/minimap2/align/main'
include { PORECHOP_ABI                           } from '../../modules/nf-core/porechop/abi/main'
include { CHOPPER as CHOPPER_LAMBDA              } from '../../modules/nf-core/chopper/main'
include { CHOPPER as CHOPPER_NOLAMBDA            } from '../../modules/nf-core/chopper/main'

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

        if (meta.platform == "ont"){
            if ( q20_percentage >= params.long_reads_ont_quality_threshold ) {
                return [ meta + [quality: "high"], reads]
            } else {
                return [ meta + [quality: "low"], reads]
            }
        } else if (meta.platform == "pb") {
            if ( q20_percentage >= params.long_reads_pacbio_quality_threshold ) {
                return [ meta + [quality: "high"], reads]
            } else {
                return [ meta + [quality: "low"], reads]
            }
        }
    }
    reads_quality_levels.view()

    /*********************************************************************/
    /* Remove ONT adapters from Nanopore sequences only                  */
    /*********************************************************************/
    reads_quality_levels.branch { meta, reads ->
        ont: meta.platform == "ont"
        pacbio: meta.platform == "pb"
    }.set {reads_platform}

    PORECHOP_ABI(
        reads_platform.ont
    )

    /***************************************************************************/
    /* Perform decontamination from lambda phage if requested                  */
    /***************************************************************************/

    PORECHOP_ABI.out.reads
        .branch { meta, _reads ->
            run_decontamination: meta.lambdaphage_reference != null
            skip_decontamination: meta.lambdaphage_reference == null
        }
        .set { lambdaphage_subdivided_reads }

    lambdaphage_subdivided_reads.run_decontamination
        .multiMap { meta, reads ->
            reads: [meta, reads]
            reference: file("${params.reference_genomes_folder}/${meta.lambdaphage_reference}.*")
        }
        .set { ch_lambda_decontamination_input }

    CHOPPER_LAMBDA(
        ch_lambda_decontamination_input.reads,
        ch_lambda_decontamination_input.reference
    )
    ch_versions = ch_versions.mix(CHOPPER_LAMBDA.out.versions)

    // chopper is still needed for start-of-sequence trimming
    CHOPPER_NOLAMBDA(
        lambdaphage_subdivided_reads.skip_decontamination,
        ""
    )
    ch_versions = ch_versions.mix(CHOPPER_NOLAMBDA.out.versions)

    // putting together pacbio and non-cleaned ont
    fastp_adapter_free_reads = CHOPPER_NOLAMBDA.out.fastq.mix(
        reads_platform.pb
    )

    lambdaphage_cleaned_reads = fastp_adapter_free_reads.mix(
        CHOPPER_LAMBDA.out.fastq
    )

    // TODO: add filter if too many reads are removed

    /***************************************************************************/
    /* Perform decontamination from human sequences if requested               */
    /***************************************************************************/
    lambdaphage_cleaned_reads
        .branch { meta, _reads ->
            run_decontamination: meta.human_reference != null
            skip_decontamination: meta.human_reference == null
        }
        .set { human_subdivided_reads }

    human_subdivided_reads.run_decontamination
        .multiMap { meta, reads_ ->
            reads: [meta, reads_]
            reference: [ [id:"human"], file("${params.reference_genomes_folder}/${meta.human_reference}", checkIfExists: true)]
        }
        .set { ch_human_decontamination_input }

    MINIMAP2_ALIGN_HUMAN(
        ch_human_decontamination_input.reads,
        ch_human_decontamination_input.reference,
        "human",
        "fastq", // out sequence extension
        true,    // bam_format = true enables filtering of the bam file with samtools inside minimap2 module
        false,   // no bam index extension needed
        false,   // no CIGAR in paf format
        true     // cigar_bam = true enables -L minimap2 option that is required for long reads alignment
                 // because of the CIGAR-related bug that BAM format has
                 // see explanation here https://github.com/lh3/minimap2?tab=readme-ov-file#working-with-65535-cigar-operations
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
        .multiMap { meta, reads_ ->
            reads: [meta, reads_]
            reference: [ [id:file(meta.contaminant_reference).baseName], file("${params.reference_genomes_folder}/${meta.contaminant_reference}", checkIfExists: true) ]
        }
        .set { ch_decontamination_input }

    MINIMAP2_ALIGN_HOST(
        ch_decontamination_input.reads,
        ch_decontamination_input.reference,
        "host",
        "fastq", // out sequence extension
        true,    // bam_format = true enables filtering of the bam file with samtools inside minimap2 module
        false,   // no bam index extension needed
        false,   // no CIGAR in paf format
        true     // cigar_bam = true enables -L minimap2 option that is required for long reads alignment
                 // because of the CIGAR-related bug that BAM format has
                 // see explanation here https://github.com/lh3/minimap2?tab=readme-ov-file#working-with-65535-cigar-operations
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
