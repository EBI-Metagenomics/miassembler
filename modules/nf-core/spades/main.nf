process SPADES {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spades:3.15.5--h95f258a_1' :
        'biocontainers/spades:3.15.5--h95f258a_1' }"

    input:
    tuple val(meta), path(illumina), path(pacbio), path(nanopore)
    path yml
    path hmm

    output:
    tuple val(meta), path('*.scaffolds.fa.gz')                     , optional:true, emit: scaffolds
    tuple val(meta), path('*.contigs.fa.gz')                       , optional:true, emit: contigs
    tuple val(meta), path('*.transcripts.fa.gz')                   , optional:true, emit: transcripts
    tuple val(meta), path('*.gene_clusters.fa.gz')                 , optional:true, emit: gene_clusters
    tuple val(meta), path('*.assembly_graph_with_scaffolds.gfa.gz'), optional:true, emit: gfa
    tuple val(meta), path('*.assembly_graph.fastg.gz')             , optional:true, emit: fastg
    tuple val(meta), path('params.txt')                            , optional:true, emit: params
    tuple val(meta), path('*.log')                                 , emit: log
    path  "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxmem = task.memory.toGiga()
    def illumina_reads = illumina ? ( meta.single_end ? "-s $illumina" : "-1 ${illumina[0]} -2 ${illumina[1]}" ) : ""
    def pacbio_reads = pacbio ? "--pacbio $pacbio" : ""
    def nanopore_reads = nanopore ? "--nanopore $nanopore" : ""
    def iontorrent_flag = meta.platform == "iontorrent" ? "--iontorrent" : ""
    def custom_hmms = hmm ? "--custom-hmms $hmm" : ""
    def reads = yml ? "--dataset $yml" : "$illumina_reads $pacbio_reads $nanopore_reads"
    def metaspades_arg = meta.assembler == "metaspades" ? "--meta" : ""
    // FIXME: figure out how to use spades retry mechanism
    //        the problem modifyng the flags or memory forces the expiration
    //        the nextflow cache hence a new working directory, 
    //        which then makes spades fail as if needs the "checkpointed files to resume itself"
    """
    spades.py \\
        $args \\
        $metaspades_arg \\
        --threads $task.cpus \\
        --memory $maxmem \\
        $custom_hmms \\
        $iontorrent_flag \\
        $reads \\
        -o ./
    mv spades.log ${prefix}.spades.log

    if [ -f scaffolds.fasta ]; then
        mv scaffolds.fasta ${prefix}.scaffolds.fa
        gzip -n ${prefix}.scaffolds.fa
    fi
    if [ -f contigs.fasta ]; then
        mv contigs.fasta ${prefix}.contigs.fa
        gzip -n ${prefix}.contigs.fa
    fi
    if [ -f transcripts.fasta ]; then
        mv transcripts.fasta ${prefix}.transcripts.fa
        gzip -n ${prefix}.transcripts.fa
    fi
    if [ -f assembly_graph_with_scaffolds.gfa ]; then
        mv assembly_graph_with_scaffolds.gfa ${prefix}.assembly_graph_with_scaffolds.gfa
        gzip -n ${prefix}.assembly_graph_with_scaffolds.gfa
    fi
    if [ -f assembly_graph.fastg ]; then
        mv assembly_graph.fastg ${prefix}.assembly_graph.fastg
        gzip -n ${prefix}.assembly_graph.fastg
    fi
    if [ -f gene_clusters.fasta ]; then
        mv gene_clusters.fasta ${prefix}.gene_clusters.fa
        gzip -n ${prefix}.gene_clusters.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
    END_VERSIONS
    """
}
