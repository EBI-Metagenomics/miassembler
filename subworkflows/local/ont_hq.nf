include { FLYE         } from '../../modules/nf-core/flye/main'

workflow ONT_HQ {
    take:
    qc_reads                   // [ val(meta), path(reads) ]

    main:

    def ch_versions = Channel.empty()

    FLYE(
        qc_reads,
        "--nano-hq"
    )
    ch_versions = ch_versions.mix(FLYE.out.versions)

    emit:
    contigs  = FLYE.out.fasta
    versions = ch_versions
}
