process FETCHTOOL_METADATA {
    tag "$study_accession"

    label 'process_single'

    container "microbiome-informatics/fetch-tool:v0.9.0"

    input:
    tuple val(meta), val(study_accession), val(reads_accession)
    path fetchtool_config

    output:
    tuple val(meta), path("download_folder/${study_accession}/*.txt"), emit: metadata
    tuple val(meta), env(library_strategy)                           , emit: lib_strategy
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def library_strategy = ''
    def private_study = params.private_study ? "--private" : ""

    """
    fetch-read-tool -d download_folder/ \\
    -p $study_accession \\
    --fix-desc-file \\
    -c $fetchtool_config \\
    -v $private_study $args

    library_strategy=\$(grep $reads_accession download_folder/$study_accession/${study_accession}.txt | cut -f 7)
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fetch-tool: \$(fetch-read-tool --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    mkdir -p download_folder/$study_accession/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fetch-tool: \$(fetch-read-tool --version)
    END_VERSIONS
    """
}
