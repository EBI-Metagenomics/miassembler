process PIGZ {
    label 'process_low'
    stageInMode 'copy'
    tag "${file_to_compress}"

    input:
    file(file_to_compress)

    output:
    path("*.gz"), emit: compressed
    path "versions.yml"                       , emit: versions

    script:
    """
    pigz ${file_to_compress}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(pigz --version 2>&1 | sed 's/pigz //g')
    END_VERSIONS
    """
}
