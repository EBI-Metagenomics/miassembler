// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SPADES {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bd1aa46fe63e43f2a7bcd7cafded3236b4830bd8ec78adec0c9140e5330d908f/data' :
        'community.wave.seqera.io/library/spades:4.3.0--ca99aa8b3ff910f5' }"

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
    tuple val(meta), path('*.warnings.log')                        , optional:true, emit: warnings
    tuple val(meta), path('*.spades.log')                          , emit: log
    tuple val("${task.process}"), val('spades'), eval("spades.py --version 2>&1 | sed -n 's/^.*SPAdes genome assembler v//p'"), topic: versions, emit: versions_spades

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxmem = task.memory.toGiga()
    def illumina_reads = illumina ? ( meta.single_end ? "-s $illumina" : "-1 ${illumina[0]} -2 ${illumina[1]}" ) : ""
    def pacbio_reads = pacbio ? "--pacbio $pacbio" : ""
    def nanopore_reads = nanopore ? "--nanopore $nanopore" : ""
    def custom_hmms = hmm ? "--custom-hmms $hmm" : ""
    def reads = yml ? "--dataset $yml" : "$illumina_reads $pacbio_reads $nanopore_reads"
    // On a Nextflow retry the task runs in a new work directory. To avoid assembling from scratch,
    // reuse the SPAdes output directory of the previous attempt and restart from its last checkpoint.
    // task.previousTrace only points to the previous attempt, so its spades_out symlink is resolved
    // to the original directory (attempt 3+ link straight to the directory that holds the checkpoints).
    // --restart-from does not accept reads or flags such as --meta/--iontorrent/--only-assembler
    // (SPAdes keeps them from the original run), so only --threads and --memory are updated.
    def previous_workdir = task.previousTrace ? task.previousTrace.get('workdir') : ''
    """
    if [ -n "${previous_workdir}" ] && [ -d "${previous_workdir}/spades_out/pipeline_state" ]; then
        echo "Restarting SPAdes from the last checkpoint of the previous attempt: \$(readlink -f ${previous_workdir}/spades_out)"
        ln -s "\$(readlink -f ${previous_workdir}/spades_out)" spades_out
        spades.py \\
            --restart-from last \\
            --threads $task.cpus \\
            --memory $maxmem \\
            -o spades_out
    else
        spades.py \\
            $args \\
            --threads $task.cpus \\
            --memory $maxmem \\
            $custom_hmms \\
            $reads \\
            -o spades_out
    fi

    # Copy (not move) the final SPAdes outputs into the task work directory, so spades_out stays
    # complete in case this attempt fails later and the next one needs to restart from it
    for output_file in spades.log warnings.log params.txt scaffolds.fasta contigs.fasta transcripts.fasta \\
        gene_clusters.fasta assembly_graph_with_scaffolds.gfa assembly_graph.fastg; do
        if [ -e "spades_out/\$output_file" ]; then
            cp "spades_out/\$output_file" ./
        fi
    done

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

    if [ -f warnings.log ]; then
        mv warnings.log ${prefix}.warnings.log
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.scaffolds.fa.gz
    echo "" | gzip > ${prefix}.contigs.fa.gz
    echo "" | gzip > ${prefix}.transcripts.fa.gz
    echo "" | gzip > ${prefix}.gene_clusters.fa.gz
    echo "" | gzip > ${prefix}.assembly.gfa.gz
    touch ${prefix}.spades.log
    touch ${prefix}.warnings.log
    """
}
