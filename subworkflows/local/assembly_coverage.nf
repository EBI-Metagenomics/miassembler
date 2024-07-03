include { BWAMEM2_INDEX                        } from '../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM as BWAMEM2_MEM_COVERAGE  } from '../../modules/ebi-metagenomics/bwamem2/mem/main'
include { SAMTOOLS_IDXSTATS                    } from '../../modules/nf-core/samtools/idxstats/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS } from '../../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'

workflow ASSEMBLY_COVERAGE {

    take:
    reads_assembly   // [ val(meta), path(reads), path(assembly_fasta) ]

    main:

    ch_versions = Channel.empty()

    reads = reads_assembly.map { meta, reads, _ -> [meta, reads]}
    assembly = reads_assembly.map { meta, _, assembly -> [meta, assembly]}

    BWAMEM2_INDEX(
        assembly
    )

    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

    bwa_index_reads = BWAMEM2_INDEX.out.index.join( reads )

    BWAMEM2_MEM_COVERAGE(
        bwa_index_reads
    )

    ch_versions = ch_versions.mix(BWAMEM2_MEM_COVERAGE.out.versions)

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS(
        BWAMEM2_MEM_COVERAGE.out.bam
    )

    ch_versions = ch_versions.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions)

    SAMTOOLS_IDXSTATS(
        BWAMEM2_MEM_COVERAGE.out.bam
    )

    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)

    emit:
    coverage_depth     = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
    samtools_idxstats  = SAMTOOLS_IDXSTATS.out.idxstats
    versions           = ch_versions
}
