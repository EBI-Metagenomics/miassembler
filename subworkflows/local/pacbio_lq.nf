include { CANU as CANU_PACBIO       } from '../../modules/nf-core/canu/main'
include { FLYE                      } from '../../modules/nf-core/flye/main'
include { MINIMAP as MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include { RACON                     } from '../../modules/nf-core/racon/main'

workflow PACBIO_LQ {
    take:
    reads                   // [ val(meta), path(reads) ]

    main:
    CANU_PACBIO(
        reads,
        "-pacbio",
        "5m"
    )
    ch_versions = ch_versions.mix(CANU_PACBIO.out.versions)

    FLYE(
        CANU_PACBIO.out.corrected_trimmed_reads,
        "--pacbio-raw"
    )
    ch_versions = ch_versions.mix(FLYE.out.versions)

    RACON(
        tuple(CANU_PACBIO.out.corrected_trimmed_reads, FLYE.out.fasta, MINIMAP.out.paf)
    )
    ch_versions = ch_versions.mix(RACON.out.versions)

    // TODO: CALL FOR BOTH HOST AND HUMAN
    MINIMAP(
        CANU_PACBIO.out.corrected_trimmed_reads,
        RACON.out.improved_assembly,
        "",
        false,
        "bai",
        false,
        false
    )
    ch_versions = ch_versions.mix(MINIMAP.out.versions)

    emit:
    contigs  = RACON.out.improved_assembly
    versions = ch_versions
}
