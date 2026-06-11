process BWAMEM2DECONTNOBAMS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bwa-mem2_samtools:0d13d0f418bae5da' :
        'community.wave.seqera.io/library/bwa-mem2_samtools:df7fb96cc09814fa' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*{_1,_2,_interleaved}.fq.gz"), emit: decont_reads
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref_prefix = task.ext.ref_prefix ?: "${meta2.id}"
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    if [[ "${meta.single_end}" == "true" ]]; then
        bwa-mem2 \\
            mem \\
            -M \\
            -t $task.cpus \\
            \$INDEX \\
            $reads \\
            | samtools view -@ ${task.cpus} -f 4 -F 256 -uS - \\
            | samtools sort -@ ${task.cpus} -n -O bam - \\
            | samtools bam2fq -@ $task.cpus - | gzip --no-name > ${ref_prefix}_${prefix}_interleaved.fq.gz
    else
        bwa-mem2 \\
            mem \\
            -M \\
            -t $task.cpus \\
            \$INDEX \\
            $reads \\
            | samtools view -@ ${task.cpus} -f 4 -F 256 -uS - \\
            | samtools sort -@ ${task.cpus} -n -O bam - \\
            | samtools bam2fq -@ ${task.cpus} -1 ${ref_prefix}_${prefix}_1.fq.gz -2 ${ref_prefix}_${prefix}_2.fq.gz -0 /dev/null -s /dev/null
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2> /dev/null)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
