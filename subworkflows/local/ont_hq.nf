include { METAMDBG_ASM } from '../../modules/nf-core/metamdbg/asm/main'

workflow ONT_HQ {
    take:
    qc_reads                   // [ val(meta), path(reads) ]

    main:

    def ch_versions = channel.empty()

    METAMDBG_ASM(
        qc_reads,
        "ont"
    )
    ch_versions = ch_versions.mix(METAMDBG_ASM.out.versions_metamdbg)

    emit:
    contigs  = METAMDBG_ASM.out.contigs
    versions = ch_versions
}
