process PUBLISH_FILE {

    publishDir(
        path: "${params.outdir}/$output_dir",
        mode: params.publish_dir_mode
    )

    input:
    path(contigs)
    val(output_dir)

    script:
    """
    echo "$output_dir"
    """
}
