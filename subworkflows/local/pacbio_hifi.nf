// placeholder for hifiadapterfilt
include { HIFIADAPTERFILT } from '../../modules/...hifiadapterfilt/main'
include { FLYE            } from '../../modules/nf-core/flye/main'

workflow PACBIO_HIFI {
    take:
    reads                   // [ val(meta), path(reads) ]

    main:
    HIFIADAPTERFILT(
        reads
    )
    ch_versions = ch_versions.mix(HIFIADAPTERFILT.out.versions)

    FLYE(
        CANU_ONT.out.corrected_trimmed_reads,
        "--pacbio-hifi"
    )
    ch_versions = ch_versions.mix(FLYE.out.versions)

    emit:
    contigs  = FLYE.out.fasta
    versions = ch_versions
}
