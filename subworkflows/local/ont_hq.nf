include { PORECHOP_ABI } from '../../modules/nf-core/porechop/abi/main'
include { FLYE         } from '../../modules/nf-core/flye/main'

workflow ONT_HQ {
    take:
    reads                   // [ val(meta), path(reads) ]

    main:

    def ch_versions = Channel.empty()

    PORECHOP_ABI(
        reads
    )
    ch_versions = ch_versions.mix(PORECHOP_ABI.out.versions)

    FLYE(
        PORECHOP_ABI.out.reads,
        "--nano-hq"
    )
    ch_versions = ch_versions.mix(FLYE.out.versions)

    emit:
    contigs  = FLYE.out.fasta
    versions = ch_versions
}
