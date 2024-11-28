include { MINIMAP2_ALIGN as MINIMAP_COVERAGE   } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_IDXSTATS                    } from '../../modules/nf-core/samtools/idxstats/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS } from '../../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'

workflow LONG_READS_ASSEMBLY_COVERAGE {

    take:
    assembly_reads   // [ val(meta), path(assembly_fasta), path(reads) ]

    main:

    ch_versions = Channel.empty()

    reads = assembly_reads.map { meta, _, reads -> [meta, reads] }
    assembly = assembly_reads.map { meta, assembly, _ -> [meta, assembly] }

    MINIMAP_COVERAGE(
        reads,
        assembly,
        "coverage",
        "bam",      // out sequence extension
        true,       // no bam format
        "bai",      // extension needed
        false,      // no CIGAR in bam format
        false       // no CIGAR in paf format
    )

    ch_versions = ch_versions.mix(MINIMAP_COVERAGE.out.versions)

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS(
        MINIMAP_COVERAGE.out.bam.join(MINIMAP_COVERAGE.out.index)
    )

    ch_versions = ch_versions.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions)

    SAMTOOLS_IDXSTATS(
        MINIMAP_COVERAGE.out.bam.join(MINIMAP_COVERAGE.out.index)
    )

    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)

    emit:
    coverage_depth     = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
    samtools_idxstats  = SAMTOOLS_IDXSTATS.out.idxstats
    versions           = ch_versions
}
