process RAW_READ_QUALITY_CHECK {
    tag "$reads_accession"
    label 'process_single'

    input:
    tuple val(meta), path(fastp_json)

    output:
    env(quality)       , emit: quality
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    quality=\$(check_raw_quality.py -j ${fastp_json})

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
