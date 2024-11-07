process FETCHTOOL_READS {
    tag "$study_accession - $reads_accession"

    label 'process_single'

    container "quay.io/microbiome-informatics/fetch-tool:v1.0.2"

    input:
    tuple val(meta), val(study_accession), val(reads_accession)
    path fetchtool_config

    output:
    tuple val(meta), path("download_folder/${study_accession}/raw/${reads_accession}*.fastq.gz"), env(library_strategy), env(library_layout), env(platform), emit: reads
    // The '_mqc.' is for multiQC
    tuple val(meta), path("download_folder/${study_accession}/${study_accession}.txt")      , emit: metadata_tsv
    path "versions.yml"                                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${reads_accession}"
    def private_study = params.private_study ? "--private" : ""
    """
    fetch-read-tool -d download_folder/ \\
    -p ${study_accession} \\
    -ru ${reads_accession} \\
    -c ${fetchtool_config} \\
    -v ${private_study} ${args}

    library_strategy=\$(echo "\$(grep ${reads_accession} download_folder/${study_accession}/${study_accession}.txt | cut -f 7)" | tr '[:upper:]' '[:lower:]')
    library_layout=\$(echo "\$(grep ${reads_accession} download_folder/${study_accession}/${study_accession}.txt | cut -f 5)" | tr '[:upper:]' '[:lower:]')

    export metadata_platform=\$(echo "\$(grep ${reads_accession} download_folder/${study_accession}/${study_accession}.txt | cut -f 8)" | tr '[:upper:]' '[:lower:]')
    if [[ \$metadata_platform == "minion" || \$metadata_platform == "promethion" || \$metadata_platform == "gridion" ]]; then
        platform="ont"
    elif [[ \$metadata_platform == "pacbio rs" || \$metadata_platform == "pacbio rs ii" ]]; then
        platform="pacbio"
    else
        platform=\$metadata_platform
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fetch-tool: \$(fetch-read-tool --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$reads_accession"
    """
    mkdir -p download_folder/$study_accession/raw/

    touch download_folder/$study_accession/raw/${reads_accession}.fastq.gz

    touch download_folder/${study_accession}/fetch_tool_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fetch-tool: \$(fetch-read-tool --version)
    END_VERSIONS
    """
}