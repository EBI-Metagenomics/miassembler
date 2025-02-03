include { LONG_READS_COVERAGE         } from '../../modules/local/long_reads_coverage/main'
include { CALCULATE_ASSEMBLY_COVERAGE } from '../../modules/local/calculate_assembly_coverage'

workflow LONG_READS_ASSEMBLY_COVERAGE {

    take:
    assembly_reads   // [ val(meta), path(assembly_fasta), path(reads) ]
    fastp_json       // [ val(meta), path(fasp_json) ]

    main:

    ch_versions = Channel.empty()

    LONG_READS_COVERAGE(
        assembly_reads
    )
    ch_versions = ch_versions.mix(LONG_READS_COVERAGE.out.versions)

    // This snippet allows to only use relevant meta values from the two, 
    // since meta of the depth file contains way more fields than the
    // previous one, preventing a direct join from working
    def fastp = fastp_json.map { meta, json_file ->
        key = meta.subMap('id', 'platform')
        return [key, json_file]
    }

    def depth = LONG_READS_COVERAGE.out.depth.map { meta, depth_file ->
        key = meta.subMap('id', 'platform')
        key2 = meta.subMap('assembler', 'assembler_version')
        return [key, key2, depth_file]
    }

    def depth_fastp_json = depth.join(fastp).map{ meta, meta2, json_file, depth_file ->
        return [meta + meta2, json_file, depth_file]
    }

    // This process calculates a single coverage and coverage depth value for the whole assembly //
    CALCULATE_ASSEMBLY_COVERAGE(
        depth_fastp_json
    )
    
    ch_versions = ch_versions.mix(CALCULATE_ASSEMBLY_COVERAGE.out.versions)

    emit:
    coverage_depth         = LONG_READS_COVERAGE.out.depth
    samtools_idxstats      = LONG_READS_COVERAGE.out.idxstats
    assembly_coverage_json = CALCULATE_ASSEMBLY_COVERAGE.out.assembly_coverage_json
    versions               = ch_versions
}
