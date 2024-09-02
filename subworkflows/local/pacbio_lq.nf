include { CANU as CANU_PACBIO                  } from '../../modules/nf-core/canu/main'

workflow PACBIO_LQ {
    take:
    reads                   // [ val(meta), path(reads) ]

    main:
    CANU_PACBIO(
        reads,
        "-pacbio",
        "5m"
    )
    CANU_PACBIO.out.corrected_reads.view()
}
