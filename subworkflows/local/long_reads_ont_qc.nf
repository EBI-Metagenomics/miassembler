include { PORECHOP_ABI                           } from '../../modules/nf-core/porechop/abi/main'
include { CHOPPER as CHOPPER_LAMBDA              } from '../../modules/nf-core/chopper/main'
include { CHOPPER as CHOPPER_NOLAMBDA            } from '../../modules/nf-core/chopper/main'

workflow LONG_READS_ONT_QC {

    take:
    ont_reads

    main:
    ch_versions = Channel.empty()

    PORECHOP_ABI(
        ont_reads
    )
    ch_versions = ch_versions.mix(PORECHOP_ABI.out.versions)

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

    lambdaphage_cleaned_reads = CHOPPER_NOLAMBDA.out.fastq.mix(
        CHOPPER_LAMBDA.out.fastq
    )

    emit:
    ont_qc_reads = lambdaphage_cleaned_reads
    versions     = ch_versions
}