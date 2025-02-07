include { SHORT_READS_INDEX_FASTA     } from '../../modules/local/short_reads_coverage/main'
include { SHORT_READS_COVERAGE        } from '../../modules/local/short_reads_coverage/main'
include { CALCULATE_ASSEMBLY_COVERAGE } from '../../modules/local/calculate_assembly_coverage'

workflow SHORT_READS_ASSEMBLY_COVERAGE {

    take:
    assembly_reads   // [ val(meta), path(assembly_fasta), path(reads) ]
    fastp_json       // [ val(meta), path(fasp_json) ]

    main:

    ch_versions = Channel.empty()

    def reads = assembly_reads.map { meta, __, reads -> [meta, reads] }
    def assembly = assembly_reads.map { meta, assembly, __ -> [meta, assembly] }

    SHORT_READS_INDEX_FASTA(
        assembly
    )

    reads_assembly_index = reads.join(SHORT_READS_INDEX_FASTA.out.fasta_with_index)

    SHORT_READS_COVERAGE(
        reads_assembly_index
    )

    ch_versions = ch_versions.mix( SHORT_READS_INDEX_FASTA.out.versions.first() )
    ch_versions = ch_versions.mix( SHORT_READS_COVERAGE.out.versions.first() )

    // This process calculates a single coverage and coverage depth value for the whole assembly //
    CALCULATE_ASSEMBLY_COVERAGE(
        SHORT_READS_COVERAGE.out.depth.join ( fastp_json )
    )

    ch_versions = ch_versions.mix(CALCULATE_ASSEMBLY_COVERAGE.out.versions)

    emit:
    coverage_depth_summary        = SHORT_READS_COVERAGE.out.depth
    samtools_idxstats             = SHORT_READS_COVERAGE.out.idxstats
    assembly_coverage_json        = CALCULATE_ASSEMBLY_COVERAGE.out.assembly_coverage_json
    versions                      = ch_versions
}
