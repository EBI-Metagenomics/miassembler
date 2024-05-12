include { ASSEMBLER_VERSION } from '../../modules/local/assembler_version'
include { QUAST             } from '../../modules/nf-core/quast/main'

workflow ASSEMBLY_STATS {

    take:
    filtered_contigs        // [ val(meta), path(contigs) ]
    assembler_log

    main:
    ch_versions = Channel.empty()

    ASSEMBLER_VERSION(
        assembler_log
    )
    
    /* The QUAST module was modified to run metaQUAST instead */
    QUAST(
        filtered_contigs,
        [ [], [] ], // reference
        [ [], [] ]  // gff
    )

    ch_versions = ch_versions.mix(QUAST.out.versions)


    /*
    Things to collect:
    - assembly stats: either from megahit or metaspades"
    - 
    */

    emit:
    quast_results = QUAST.out.results
    versions = ch_versions
}
