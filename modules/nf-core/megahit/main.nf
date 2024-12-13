process MEGAHIT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0' :
        'biocontainers/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("megahit_out/*.contigs.fa.gz")                            , emit: contigs
    tuple val(meta), path("megahit_out/*.fastg.gz")                                 , emit: graph
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
    if (meta.single_end) {
        """
        megahit \\
            -r ${reads} \\
            -t $task.cpus \\
            $args \\
            --out-prefix $prefix

        if [ ! -s megahit_out/*.fa ]; then
            echo "No contigs assembled" | tee /dev/stderr
            exit 1
        fi

        export kmax=\$(grep "k-max reset to" megahit_out/*log | rev | cut -f 2 -d ' ' | rev)
        megahit_toolkit contig2fastg \$kmax megahit_out/*.fa > megahit_out/${prefix}.contigs.fastg

        pigz \\
            --no-name \\
            -p $task.cpus \\
            $args2 \\
            megahit_out/*.fa \\
            megahit_out/intermediate_contigs/*.fa \\
            megahit_out/*.fastg

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
            --out-prefix $prefix

        if [ ! -s megahit_out/*.fa ]; then
            echo "No contigs assembled" | tee /dev/stderr
            exit 1
        fi

        export kmax=\$(grep "k-max reset to" megahit_out/*log | rev | cut -f 2 -d ' ' | rev)
        megahit_toolkit contig2fastg \$kmax megahit_out/*.fa > megahit_out/${prefix}.contigs.fastg

        pigz \\
            --no-name \\
            -p $task.cpus \\
            $args2 \\
            megahit_out/*.fa \\
            megahit_out/intermediate_contigs/*.fa \\
            megahit_out/*.fastg

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
        END_VERSIONS
        """
    }
}
