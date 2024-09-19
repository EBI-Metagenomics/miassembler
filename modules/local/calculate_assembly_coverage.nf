process CALCULATE_ASSEMBLY_COVERAGE {

    tag "$meta.id"

    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'quay.io/biocontainers/python:3.10' }"

    input:
    tuple val(meta), path(jgi_summary_tsv_gz), path(fastp_json)

    output:
    tuple val(meta), path("${meta.id}_coverage.json"), emit: assembly_coverage_json
    val("versions.yml"),                                        emit: version

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    calculate_assembly_covarege.py -j ${jgi_summary_tsv_gz} -f ${fastp_json} -o ${prefix}_coverage.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_coverage_data.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
