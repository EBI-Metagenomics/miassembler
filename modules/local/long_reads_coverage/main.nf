process LONG_READS_COVERAGE {

    /*
    This module aligns reads to the reference and produces output in the form of FASTQ.GZ, BAM and BAM.BAI files.
    It also runs depth generation using jgi_summarize_bam_contig_depths
    */

    label 'process_high'

    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/21/21ad589eae1b0fc318f4c086886799c96370adfb6dc8561cb62265fe527d9fb2/data':
        'community.wave.seqera.io/library/metabat2_minimap2_samtools:befe455c29d07c61' }"

    input:
    tuple val(meta), path(assembly_fasta), path(reads)

    output:
    tuple val(meta), path("*.tsv.gz")                     , emit: depth
    tuple val(meta), path("*.idxstats")                   , emit: idxstats
    path "versions.yml"                                   , emit: versions

    script:

    def prefix = task.ext.prefix ?: "${meta.id}"

    def jgi_summarize_bam_contig_depths_args = task.ext.jgi_summarize_bam_contig_depths_args ?: ''

    """
    mkdir -p output
    echo " ---> mapping files to assembly"

    minimap2 \\
        -t ${task.cpus} \\
        -x map-${meta.platform} \\
        ${assembly_fasta} \\
        ${reads} \\
        -a | samtools sort -@ ${task.cpus-1} \\
        -O bam - -o output/${meta.id}_sorted.bam

    echo " ---> samtools index sorted bam"
    samtools index -@ ${task.cpus} output/${meta.id}_sorted.bam

    echo " ---> samtools idxstats sorted bam"
    samtools idxstats --threads ${task.cpus-1} output/${meta.id}_sorted.bam > ${meta.id}.assembly.idxstats

    echo " ---> depth generation"
    jgi_summarize_bam_contig_depths \
        --outputDepth ${prefix}_coverage_depth_summary.tsv \
        $jgi_summarize_bam_contig_depths_args \
        output/${meta.id}_sorted.bam
    bgzip --threads $task.cpus ${prefix}_coverage_depth_summary.tsv

    rm -rf fasta_outdir output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2> /dev/null)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        metabat2: \$( metabat2 --help 2>&1 | head -n 2 | tail -n 1 | sed -n 's/.*version \\([^;]*\\);.*/\\1/p' )
    END_VERSIONS
    """
}
