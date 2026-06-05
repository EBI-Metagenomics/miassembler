include { PROOVFRAME_FIX } from '../../modules/nf-core/proovframe/fix/main'
include { PROOVFRAME_MAP } from '../../modules/nf-core/proovframe/map/main'

workflow FRAMESHIFT_CORRECTION {
    take:
    contigs                   // [ val(meta), path(contigs) ]

    main:
    def ch_versions = channel.empty()

    PROOVFRAME_MAP(
        contigs,
        [ [],[] ],
        [ [id:'diamond_db'], "${params.diamond_db}" ]
    )
    ch_versions = ch_versions.mix(PROOVFRAME_MAP.out.versions)

    pf_fix_input = contigs
        .join(PROOVFRAME_MAP.out.tsv)
        .multiMap { meta, assembly_contigs, tsv ->
            contigs:    [meta, assembly_contigs]
            pf_map_out: [meta, tsv]
    }

    PROOVFRAME_FIX(
        pf_fix_input.contigs,
        pf_fix_input.pf_map_out
    )
    ch_versions = ch_versions.mix(PROOVFRAME_FIX.out.versions)

    emit:
    corrected_contigs  = PROOVFRAME_FIX.out.out_fa
    versions           = ch_versions
}
