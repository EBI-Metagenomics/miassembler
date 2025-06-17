include { CANU as CANU_ONT } from '../../modules/nf-core/canu/main'
include { FLYE             } from '../../modules/nf-core/flye/main'
include { MINIMAP2_ALIGN   } from '../../modules/nf-core/minimap2/align/main'
include { RACON            } from '../../modules/nf-core/racon/main'
include { MEDAKA           } from '../../modules/nf-core/medaka/main'

workflow ONT_LQ {
    take:
    reads                   // [ val(meta), path(reads) ]

    main:

    ch_versions = Channel.empty()

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
        "",     // no prefix needed for paf generation
        "",     // no extension needed for paf generation
        false,  // no bam format
        false,  // no bam index extension needed for paf generation
        false,  // no CIGAR in paf format
        false   // no CIGAR in bam format
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    reads_flye_contigs_paf = CANU_ONT.out.corrected_trimmed_reads
            .join(FLYE.out.fasta)
            .join(MINIMAP2_ALIGN.out.paf)

    RACON(
        reads_flye_contigs_paf
    )
    ch_versions = ch_versions.mix(RACON.out.versions)

    reads_racon_contigs = CANU_ONT.out.corrected_trimmed_reads.join(RACON.out.improved_assembly)

    MEDAKA(
        reads_racon_contigs
    )
    ch_versions = ch_versions.mix(MEDAKA.out.versions)

    emit:
    contigs  = MEDAKA.out.assembly
    versions = ch_versions
}
