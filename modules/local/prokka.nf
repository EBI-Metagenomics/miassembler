process PROKKA {

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl526_0' :
        'biocontainers/prokka:1.14.6--pl526_0' }"

    input:
    tuple val(meta), path(assembly_file)

    output:
    tuple val(meta), path('prokka_out/contigs.gbk'), emit: prokka_gbk
	tuple val(meta), path('prokka_out/contigs.gff'), emit: prokka_gff
	tuple val(meta), path('prokka_out/contigs.faa'), emit: prokka_faa
	tuple val(meta), path('prokka_out/contigs.fna'), emit: prokka_fna
    path "versions.yml",                             emit: versions

    script:
    def is_compressed = assembly_file.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? assembly_file.getBaseName() : assembly_file

    if(assembly_file.size() > 0)
        """
        if [ "${is_compressed}" == "true" ]; then
            gzip -c -d ${assembly_file} > ${fasta_name}
        fi

        # TMP folder issues in Prokka - https://github.com/tseemann/prokka/issues/402
        export TMPDIR="\$PWD/tmp"
        mkdir -p "\$PWD/tmp"
        # Disable the Java VM performane gathering tool, for improved performance
        export JAVA_TOOL_OPTIONS="-XX:-UsePerfData"

        prokka --outdir prokka_out \
        --prefix contigs \
        --cpus ${task.cpus} \
        --metagenome \
        ${fasta_name}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            prokka: \$(prokka --version)
        END_VERSIONS
        """
    else
        """
        echo 'PROKKA dir empty due to empty input... generating dummy files'
        mkdir prokka_out
        touch prokka_out/contigs.gbk
        """
}