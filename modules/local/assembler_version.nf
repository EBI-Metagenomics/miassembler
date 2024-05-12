process ASSEMBLER_VERSION {
    input:
    tuple val(meta), val(log_file)
    
    output:
    tuple val(meta), env(assembler_sw)                        , emit: assembler
    tuple val(meta), path("${meta.id}/assembler_version.json"), emit: stats_file
    path "versions.yml"                                       , emit: versions
    
    script:
    def read_accession = "${meta.id}"
    def assembler_sw = ''
    def assembler_vs = ''

    """
    mkdir $read_accession
    
    get_assembler_version.py ${log_file} > ${read_accession}/assembler_version.json
    
    assembler_sw=\$(head -n 1 ${read_accession}/assembler_version.json | cut -f 2 -d ':' | cut -f 1 -d ',')
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
