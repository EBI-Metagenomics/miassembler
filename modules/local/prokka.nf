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
    if(assembly_file.size() > 0)
        """
        # TMP folder issues in Prokka - https://github.com/tseemann/prokka/issues/402
        export TMPDIR="\$PWD/tmp"
        mkdir -p "\$PWD/tmp"
        # Disable the Java VM performane gathering tool, for improved performance
        export JAVA_TOOL_OPTIONS="-XX:-UsePerfData"

        gunzip -c ${assembly_file} > unzipped_assembly.fasta

        prokka --outdir prokka_out \
        --prefix contigs \
        --cpus ${task.cpus} \
        --metagenome \
        unzipped_assembly.fasta

        rm unzipped_assembly.fasta

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