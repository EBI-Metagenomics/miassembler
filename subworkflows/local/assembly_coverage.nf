include { BWAMEM2_INDEX                        } from '../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM                          } from '../../modules/ebi-metagenomics/bwamem2/mem/main'
include { SAMTOOLS_IDXSTATS                    } from '../../modules/nf-core/samtools/idxstats/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS } from '../../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'

workflow ASSEMBLY_COVERAGE {

    take:
    reads            // [ val(meta), path(reads) ]
    assembly         // [ val(meta), path(assembly_fasta) ]

    main:

    ch_versions = Channel.empty()

    BWAMEM2_INDEX(
        assembly
    )

    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

    BWAMEM2_MEM(
        reads,
        BWAMEM2_INDEX.out.index
    )
   
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS(
        BWAMEM2_MEM.out.bam
    )

    ch_versions = ch_versions.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions)

    SAMTOOLS_IDXSTATS(
        BWAMEM2_MEM.out.bam
    )

    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)

    emit:
    coverage_depth     = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
    samtools_idxstats  = SAMTOOLS_IDXSTATS.out.idxstats
    versions           = ch_versions
}
