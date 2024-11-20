include { CANU as CANU_ONT } from '../../modules/nf-core/canu/main'
include { FLYE             } from '../../modules/nf-core/flye/main'
include { MINIMAP2_ALIGN   } from '../../modules/nf-core/minimap2/align/main'
include { RACON            } from '../../modules/nf-core/racon/main'
include { MEDAKA           } from '../../modules/nf-core/medaka/main'

workflow ONT_LQ {
    take:
    reads                   // [ val(meta), path(reads) ]

    main:
    CANU_ONT(
        reads,
        "-nanopore",
        "5m"
    )
    ch_versions = ch_versions.mix(CANU_ONT.out.versions)

    FLYE(
        CANU_ONT.out.corrected_trimmed_reads,
        "--nano-raw"
    )
    ch_versions = ch_versions.mix(FLYE.out.versions)

    MINIMAP2_ALIGN(
        CANU_ONT.out.corrected_trimmed_reads,
        FLYE.out.fasta,
        "",
        false,
        bai,
        false,
        false
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    RACON(
        tuple(CANU_ONT.out.corrected_trimmed_reads,
              FLYE.out.fasta,
              MINIMAP2_ALIGN.out.paf)
    )
    ch_versions = ch_versions.mix(RACON.out.versions)

    MEDAKA(
        tuple(CANU_ONT.out.corrected_trimmed_reads,
        RACON.out.improved_assembly)
    )
    ch_versions = ch_versions.mix(MEDAKA.out.versions)

    emit:
    contigs  = MEDAKA.out.improved_assembly
    versions = ch_versions
}
