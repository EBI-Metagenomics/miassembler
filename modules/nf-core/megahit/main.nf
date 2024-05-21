process MEGAHIT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0' :
        'biocontainers/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0' }"

    publishDir(
        path: "${params.outdir}",
        mode: params.publish_dir_mode,
        failOnError: true,
        pattern: "megahit_out/*.fa*.gz",
        saveAs: {
            filename -> {
                def output_file = new File(filename);
                def studyAccessionPrefix = params.study_accession.substring(0, 7);
                def readsAccessionPrefix = params.reads_accession.substring(0, 7);
                return "${studyAccessionPrefix}/${params.study_accession}/${readsAccessionPrefix}/${params.reads_accession}/assembly/${meta.assembler}/${meta.assembler_version}/${output_file.name}";
                }
        }
    )

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("megahit_out/*.contigs.fa.gz")                            , emit: contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.contigs.fa.gz")      , emit: k_contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.addi.fa.gz")         , emit: addi_contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.local.fa.gz")        , emit: local_contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.final.contigs.fa.gz"), emit: kfinal_contigs
    path "versions.yml"                                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def restart = ""
    if (task.attempt > 1) {
        // Set of extra flags to restart the assembly process
        restart = "--continue"
    }

    if (meta.single_end) {
        """
        megahit \\
            -r ${reads} \\
            -t $task.cpus \\
            $args \\
            $restart \\
            --out-prefix $prefix

        pigz \\
            --no-name \\
            -p $task.cpus \\
            $args2 \\
            megahit_out/*.fa \\
            megahit_out/intermediate_contigs/*.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
        END_VERSIONS
        """
    } else {
        """
        megahit \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -t $task.cpus \\
            $args \\
            $restart \\
            --out-prefix $prefix

        pigz \\
            --no-name \\
            -p $task.cpus \\
            $args2 \\
            megahit_out/*.fa \\
            megahit_out/intermediate_contigs/*.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
        END_VERSIONS
        """
    }
}
