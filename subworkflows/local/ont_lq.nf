include { CANU } from '../../modules/nf-core/canu/main'

workflow ONT_LQ {
    take:
    reads                   // [ val(meta), path(reads) ]

    main:
    CANU_ONT(
        reads,
        "-nanopore",
        "5m"
    )
    CANU_ONT.out.corrected_trimmed_reads.view()

    // temporary just to test the module
    emit:
    contigs = CANU_ONT.out.corrected_trimmed_reads
}
