process RENAME_MEGAHIT_HEADERS {

    label 'process_low'
    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("*.f*.gz"),        emit: contigs_renamed_headers
    path "versions.yml",                     emit: versions

    script:
    """
    rename_megahit_headers.py --input $contigs --output "${meta.id}.renamed.fasta"
    gzip ${meta.id}.renamed.fasta
    mv ${meta.id}.renamed.fasta.gz ${contigs.baseName}
    rm ${meta.id}.renamed.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
