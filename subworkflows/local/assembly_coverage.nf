include { BWAMEM2_INDEX                        } from '../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM                          } from '../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_INDEX                       } from '../modules/nf-core/samtools/index/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS } from '../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'

workflow ASSEMBLY_COVERAGE {

    take:
    reads            // [ val(meta), path(reads) ]
    assembly         // [ val(meta), path(assembly_fasta) ]

    main:

    ch_versions = Channel.empty()

    BWAMEM2_INDEX(
        assembly
    )

    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())

    BWAMEM2_MEM(
        assembly.join( reads ),
        BWAMEM2_INDEX.out.index,
        true // sort BAM
    )

    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    SAMTOOLS_INDEX(
        BWAMEM2_MEM.out.bam
    )

    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS(
        assembly.join( BWAMEM2_MEM.out.bam ).join( SAMTOOLS_INDEX.out.bai )
    )

    ch_versions = ch_versions.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions.first())

    output:
    coverage_depth   = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
    versions         = ch_versions
}
