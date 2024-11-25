process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)
    val prefix2
    val fa_fq_format
    val cigar_paf_format

    output:
    tuple val(meta), path("*.minimap*")                  , optional: true, emit: filtered_output
    tuple val(meta), path("*.paf")                       , optional: true, emit: paf
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def map_mode = "${meta.platform}" ? "-x map-${meta.platform}" : ''
    def output = fa_fq_format ? "-a | samtools ${fa_fq_format} -f 4 | gzip > ${prefix}.${prefix2}.minimap.${fa_fq_format}.gz" : "-o ${prefix}.paf"
    def cigar_paf = cigar_paf_format ? "-c" : ''
    def query = reads
    def target = reference

    """
    minimap2 \\
        $args \\
        -t $task.cpus \\
        $map_mode \\
        $target \\
        $query \\
        $cigar_paf \\
        $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_file = fa_fq_format ? "${prefix}.${fa_fq_format}" : "${prefix}.paf"

    """
    touch $output_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
