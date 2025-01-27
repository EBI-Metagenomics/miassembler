process HIFIADAPTERFILT {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/microbiome-informatics/hifiadapterfilt:v1.0.1'

    containerOptions {
        if (workflow.containerEngine == 'singularity') {
            return "--bind ${moduleDir}/tests/fixtures/:/data/"
        } else {
            return "-v ${moduleDir}/tests/fixtures/:/data/"
        }
    }

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.filt.fastq.gz"), emit: filt
    path("*.contaminant.blastout")          , emit: blast_search
    path("*.stats")                         , emit: stats
    path("*.blocklist")                     , emit: headers
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    hifiadapterfilt.sh ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hifiadapterfilt: \$(hifiadapterfilt.sh -v)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.filt.fastq.gz
    touch ${prefix}.contaminant.blastout
    touch ${prefix}.stats
    touch ${prefix}.blocklist

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hifiadapterfilt: \$(hifiadapterfilt.sh -v)
    END_VERSIONS
    """
}
