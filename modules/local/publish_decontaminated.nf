process PUBLISH_DECONTAMINATED {

    input:
    tuple val(meta), path(decontaminated_contigs)

    output:
    path("decontaminated.txt")

    script:
    """
    mv ${decontaminated_contigs} decontaminated.txt
    """
}
