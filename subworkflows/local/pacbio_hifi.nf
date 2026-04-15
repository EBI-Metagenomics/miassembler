include { HIFIADAPTERFILT } from '../../modules/local/hifiadapterfilt/main'
include { METAMDBG_ASM    } from '../../modules/nf-core/metamdbg/asm/main'

workflow PACBIO_HIFI {
    take:
    reads                   // [ val(meta), path(reads) ]

    main:

    ch_versions = channel.empty()

    HIFIADAPTERFILT(
        reads
    )
    ch_versions = ch_versions.mix(HIFIADAPTERFILT.out.versions)

    METAMDBG_ASM(
        HIFIADAPTERFILT.out.filt,
        "hifi"
    )
    ch_versions = ch_versions.mix(METAMDBG_ASM.out.versions)

    emit:
    contigs  = METAMDBG_ASM.out.contigs
    versions = ch_versions
}
