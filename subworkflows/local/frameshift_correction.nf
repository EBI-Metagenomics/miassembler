include { PROOVFRAME } from '../../modules/nf-core/proovframe/main'

workflow FRAMESHIFT_CORRECTION {
    take:
    contigs                   // [ val(meta), path(contigs) ]

    // TODO: temporary
    main:
    PROOVFRAME(
        reads
    )
    ch_versions = ch_versions.mix(PROOVFRAME.out.versions)

    emit:
    corrected_contigs  = PROOVFRAME.out.fasta
    versions = ch_versions
}
