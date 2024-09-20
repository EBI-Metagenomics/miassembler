include { CALCULATE_ASSEMBLY_COVERAGE          } from '../../modules/local/calculate_assembly_coverage'
include { BWAMEM2_MEM as BWAMEM2_MEM_COVERAGE  } from '../../modules/ebi-metagenomics/bwamem2/mem/main'
include { BWAMEM2_INDEX                        } from '../../modules/nf-core/bwamem2/index/main'
include { SAMTOOLS_IDXSTATS                    } from '../../modules/nf-core/samtools/idxstats/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS } from '../../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'

workflow ASSEMBLY_COVERAGE {

    take:
    assembly_reads   // [ val(meta), path(assembly_fasta), path(reads) ]
    fastp_json       // [ val(meta), path(fasp_json) ]

    main:

    ch_versions = Channel.empty()

    reads = assembly_reads.map { meta, _, reads -> [meta, reads] }
    assembly = assembly_reads.map { meta, assembly, _ -> [meta, assembly] }

    BWAMEM2_INDEX(
        assembly
    )

    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

    bwa_reads_index = reads.join( BWAMEM2_INDEX.out.index )

    BWAMEM2_MEM_COVERAGE(
        bwa_reads_index
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

    // This process calculates a single coverage and coverage depth value for the whole assembly //
    CALCULATE_ASSEMBLY_COVERAGE(
        METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth.join ( fastp_json )
    )

    ch_versions = ch_versions.mix(CALCULATE_ASSEMBLY_COVERAGE.out.versions)

    emit:
    coverage_depth_summary        = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
    samtools_idxstats             = SAMTOOLS_IDXSTATS.out.idxstats
    assembly_coverage_json        = CALCULATE_ASSEMBLY_COVERAGE.out.assembly_coverage_json
    versions                      = ch_versions
}
