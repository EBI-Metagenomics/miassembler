process SHORT_READS_INDEX_FASTA {

    label 'process_medium'
    tag "${meta.id} index ${fasta}"

    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta), path("${fasta.baseName}*.?*"), emit: fasta_with_index
    path("versions.yml"),                                        emit: versions

    script:
    """
    echo "index ref genome"
    bwa-mem2 index ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2> /dev/null)
    END_VERSIONS
    """
}

process SHORT_READS_COVERAGE {

    /*
    This module aligns reads to the reference using specified arguments and produces output in the form of FASTQ.GZ, BAM and BAM.BAI files.
    It also runs depth generation using jgi_summarize_bam_contig_depths
    */

    label 'process_medium'

    tag "${meta.id} align to ${ref_fasta}"

    container 'quay.io/microbiome-informatics/bwa_metabat:2.2.1_2.15'

    input:
    tuple val(meta), path(reads), path(ref_fasta), path(ref_fasta_index)
    val(get_depth)
    val(get_concoct_table)
    val(bed)

    output:
    tuple val(meta), path("*.txt.gz")                     , emit: depth
    tuple val(meta), path("*.tsv"), path("concoct*.fasta"), emit: concoct_data
    tuple val(meta), path("*.idxstats")                   , emit: idxstats
    path "versions.yml"                                   , emit: versions

    script:

    def prefix = task.ext.prefix ?: "${meta.id}"

    def samtools_args = task.ext.alignment_args

    def jgi_summarize_bam_contig_depths_args = task.ext.jgi_summarize_bam_contig_depths_args ?: ''

    def concoct_cut_up_fasta_args = task.ext.concoct_cut_up_fasta_args ?: ''
    def concoct_bedfile    = bed ? "-b ${prefix}.bed" : ""
    def concoct_prefix = task.ext.concoct_prefix ?: "concoct_${meta.id}"
    def concoct_coverage_table_args       = task.ext.concoct_coverage_table_args ?: ''

    """
    mkdir -p output
    echo " ---> mapping files to assembly"
    bwa-mem2 mem -M \
      -t ${task.cpus} \
      ${ref_fasta} \
      ${reads} | \
    samtools view -@ ${task.cpus} ${samtools_args} - | \
    samtools sort -@ ${task.cpus} -O bam - -o output/${meta.id}_sorted.bam

    echo " ---> samtools index sorted bam"
    samtools index -@ ${task.cpus} output/${meta.id}_sorted.bam

    echo " ---> samtools idxstats sorted bam"
    samtools idxstats --threads ${task.cpus-1} output/${meta.id}_sorted.bam > ${prefix}.assembly.idxstats

    if [[ "$get_depth" == "true" ]]; then
        echo " ---> depth generation"
        jgi_summarize_bam_contig_depths \
            --outputDepth ${prefix}.txt \
            $jgi_summarize_bam_contig_depths_args \
            output/${meta.id}_sorted.bam
        bgzip --threads $task.cpus ${prefix}.txt
    else
        touch ${prefix}.txt.gz
    fi

    rm -rf fasta_outdir output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2> /dev/null)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        metabat2: \$( metabat2 --help 2>&1 | head -n 2 | tail -n 1 | sed -n 's/.*version \\([^;]*\\);.*/\\1/p' )
        concoct: \$(echo \$(concoct --version 2>&1) | sed 's/concoct //g' )
    END_VERSIONS
    """
}
