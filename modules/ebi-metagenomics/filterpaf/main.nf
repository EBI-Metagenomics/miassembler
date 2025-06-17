process FILTERPAF {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(paf_file)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: mapped_contigs_txt
    path "versions.yml"                   , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Filter PAF by query coverage and percentage identity
    awk '
        {
            query_len = \$2;
            query_start = \$3;
            query_end = \$4;
            target = \$6;
            matching_bases = \$10;
            total_bases = \$11;

            if (target == "*") {
                next;
            }

            aligned_len = query_end - query_start;
            query_coverage = aligned_len / query_len;
            identity = matching_bases / total_bases;

            if (query_coverage >= ${params.min_qcov} && identity >= ${params.min_pid}) {
                print \$1;
            }
        }
    ' ${paf_file} > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | grep -oE '[0-9]{8}')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | grep -oE '[0-9]{8}')
    END_VERSIONS
    """
}
