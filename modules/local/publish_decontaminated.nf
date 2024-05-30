process PUBLISH_DECONTAMINATED {

    input:
    tuple val(meta), path(decontaminated_contigs)

    output:
    path("${meta.id}.txt")

    script:
    """
    mv ${decontaminated_contigs} ${meta.id}.txt
    """
}
