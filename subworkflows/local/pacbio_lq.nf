include { CANU as CANU_PACBIO } from '../../modules/nf-core/canu/main'
include { FLYE                } from '../../modules/nf-core/flye/main'
include { MINIMAP2_ALIGN      } from '../../modules/local/minimap2/align/main'
include { RACON               } from '../../modules/nf-core/racon/main'

workflow PACBIO_LQ {
    take:
    reads                   // [ val(meta), path(reads) ]

    main:

    ch_versions = Channel.empty()

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

    MINIMAP2_ALIGN(
        CANU_PACBIO.out.corrected_trimmed_reads,
        FLYE.out.fasta,
        "to_assembly",
        "",     // no extension needed for paf generation
        false,  // no bam format
        "bai",  // extension needed
        false,  // no CIGAR in bam format
        false   // no CIGAR in paf format
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    reads_flye_contigs_paf = CANU_PACBIO.out.corrected_trimmed_reads
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
