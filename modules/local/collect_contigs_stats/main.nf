process COLLECT_CONTIGS_STATS {

    label 'process_low'
    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(contigs), path(coverage), path(fastp_json)

    output:
    tuple val(meta), path("*.json"),        emit: contigs_stats
    path "versions.yml",                    emit: versions

    script:
    """
    contigs_stats.py -f $contigs -c $coverage -p $fastp_json -o "${meta.id}_contig_stats.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
