include { INDEX_FASTA        } from '../../modules/local/short_reads_coverage/main'
include { FEATURED_ALIGNMENT } from '../../modules/local/short_reads_coverage/main'

workflow SHORT_READS_ASSEMBLY_COVERAGE {

    take:
    assembly_reads   // [ val(meta), path(assembly_fasta), path(reads) ]
    fastp_json       // [ val(meta), path(fasp_json) ]

    main:

    ch_versions = Channel.empty()

    def reads = assembly_reads.map { meta, __, reads -> [meta, reads] }
    def assembly = assembly_reads.map { meta, assembly, __ -> [meta, assembly] }

    INDEX_FASTA(
        assembly
    )

    reads_assembly_index = assembly_and_reads \
        .map { meta, assembly, reads -> [ meta, reads ] } \
        .join( INDEX_FASTA.out.fasta_with_index )

    FEATURED_ALIGNMENT(
        reads_assembly_index
}

    ch_versions = ch_versions.mix( INDEX_FASTA.out.versions.first() )
    ch_versions = ch_versions.mix( FEATURED_ALIGNMENT.out.versions.first() )

    // This process calculates a single coverage and coverage depth value for the whole assembly //
    CALCULATE_ASSEMBLY_COVERAGE(
        FEATURED_ALIGNMENT.out.depth.join ( fastp_json )
    )

    ch_versions = ch_versions.mix(CALCULATE_ASSEMBLY_COVERAGE.out.versions)

    emit:
    coverage_depth_summary        = FEATURED_ALIGNMENT.out.depth
    samtools_idxstats             = FEATURED_ALIGNMENT.out.idxstats
    assembly_coverage_json        = CALCULATE_ASSEMBLY_COVERAGE.out.assembly_coverage_json
    versions                      = ch_versions
}
