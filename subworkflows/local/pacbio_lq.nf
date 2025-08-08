include { FLYE                } from '../../modules/nf-core/flye/main'
include { MINIMAP2_ALIGN      } from '../../modules/nf-core/minimap2/align/main'
include { RACON               } from '../../modules/nf-core/racon/main'

workflow PACBIO_LQ {
    take:
    qc_reads                   // [ val(meta), path(reads) ]

    main:

    ch_versions = Channel.empty()

    FLYE(
        qc_reads,
        "--pacbio-raw"
    )
    ch_versions = ch_versions.mix(FLYE.out.versions)

    MINIMAP2_ALIGN(
        qc_reads,
        FLYE.out.fasta,
        "",     // no prefix needed for paf generation
        "",     // no extension needed for paf generation
        false,  // no bam format
        false,  // no bam index extension needed for paf generation
        false,  // no CIGAR in paf format
        false   // no CIGAR in bam format
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    reads_flye_contigs_paf = qc_reads
            .join(FLYE.out.fasta)
            .join(MINIMAP2_ALIGN.out.paf)

    RACON(
        reads_flye_contigs_paf
    )
    ch_versions = ch_versions.mix(RACON.out.versions)

    emit:
    contigs  = RACON.out.improved_assembly
    versions = ch_versions
}
