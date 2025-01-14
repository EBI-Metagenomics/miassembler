include { PROOVFRAME_FIX        } from '../../modules/nf-core/proovframe/fix/main'
include { PROOVFRAME_MAP        } from '../../modules/nf-core/proovframe/map/main'

workflow FRAMESHIFT_CORRECTION {
    take:
    contigs                   // [ val(meta), path(contigs) ]

    main:
    def ch_versions = Channel.empty()

    PROOVFRAME_MAP(
        [[ id:'diamond_db' ], file(params.diamond_db, checkIfExists: true) ],
        contigs
    )
    ch_versions = ch_versions.mix(PROOVFRAME_MAP.out.versions)

    PROOVFRAME_FIX(
        contigs, 
        PROOVFRAME_MAP.out.tsv
    )
    ch_versions = ch_versions.mix(PROOVFRAME_FIX.out.versions)

    emit:
    corrected_contigs  = PROOVFRAME_FIX.out.out_fa
    versions           = ch_versions
}
