include { COLLECT_CONTIGS_STATS  } from '../../modules/local/collect_contigs_stats'
include { QUAST                  } from '../../modules/nf-core/quast/main'

workflow ASSEMBLY_STATS {

    take:
    input_data        // [ val(meta), path(contigs), path(coverage), path(fastp_json) ]

    main:
    ch_versions = Channel.empty()
    input_data.multiMap { meta, contigs, coverage, fastp_json ->
            contigs_stats: [ meta, contigs, coverage, fastp_json ]
            filtered_contigs: [ meta, contigs ]
        }.set {
            input
        }
    /*
    Stats:
        - input_read_count
        - limited_1000
        - limited_10000
        - limited_50000
        - num_contigs
        - assembly_length
        - largest_contig
        - n50
        - l50
        - coverage (value)
        - coverage_depth
    */
    COLLECT_CONTIGS_STATS( input.contigs_stats )
    ch_versions = ch_versions.mix(COLLECT_CONTIGS_STATS.out.versions)

    /* The QUAST module was modified to run metaQUAST instead */
    QUAST(
        input.filtered_contigs,
        [ [], [] ], // reference
        [ [], [] ]  // gff
    )

    ch_versions = ch_versions.mix(QUAST.out.versions)


    emit:
    quast_results = QUAST.out.results
    versions = ch_versions
}
