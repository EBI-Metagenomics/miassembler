include { PORECHOP_ABI } from '../../modules/nf-core/porechop/abi/main'

workflow ONT_HQ {
    take:
    reads                   // [ val(meta), path(reads) ]

    main:
    PORECHOP_ABI(
        reads
    )
    PORECHOP_ABI.out.reads.view()

    // temporary just to test the module
    emit:
    contigs = PORECHOP_ABI.out.reads
}
