include { HIFIADAPTERFILT } from '../../modules/ebi-metagenomics/hifiadapterfilt/main'
include { FLYE            } from '../../modules/nf-core/flye/main'

workflow PACBIO_HIFI {
    take:
    reads                   // [ val(meta), path(reads) ]

    main:

    ch_versions = Channel.empty()

    HIFIADAPTERFILT(
        reads
    )
    ch_versions = ch_versions.mix(HIFIADAPTERFILT.out.versions)

    FLYE(
        HIFIADAPTERFILT.out.filt,
        "--pacbio-hifi"
    )
    ch_versions = ch_versions.mix(FLYE.out.versions)

    emit:
    contigs  = FLYE.out.fasta
    versions = ch_versions
}
