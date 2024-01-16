process QUAST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1' :
        'biocontainers/quast:5.2.0--py39pl5321h2add14b_1' }"

    input:
    tuple val(meta) , path(consensus)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gff)

    output:
    tuple val(meta), path("${prefix}")                   , emit: results
    tuple val(meta), path("${prefix}.tsv")               , emit: tsv
    tuple val(meta), path("${prefix}_transcriptome.tsv") , optional: true , emit: transcriptome
    tuple val(meta), path("${prefix}_misassemblies.tsv") , optional: true , emit: misassemblies
    tuple val(meta), path("${prefix}_unaligned.tsv")     , optional: true , emit: unaligned
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    def features  = gff             ?  "--features $gff" : ''
    def reference = fasta           ?  "-r $fasta"       : ''
    """
    metaquast.py \\
        --output-dir $prefix \\
        $reference \\
        $features \\
        --threads $task.cpus \\
        $args \\
        ${consensus.join(' ')}

    ln -s ${prefix}/combined_reference/report.tsv ${prefix}.tsv
    [ -f  ${prefix}/combined_reference/contigs_reports/all_alignments_transcriptome.tsv ] && ln -s ${prefix}/combined_reference/contigs_reports/all_alignments_transcriptome.tsv ${prefix}_transcriptome.tsv
    [ -f  ${prefix}/combined_reference/contigs_reports/misassemblies_report.tsv         ] && ln -s ${prefix}/combined_reference/contigs_reports/misassemblies_report.tsv ${prefix}_misassemblies.tsv
    [ -f  ${prefix}/combined_reference/contigs_reports/unaligned_report.tsv             ] && ln -s ${prefix}/combined_reference/contigs_reports/unaligned_report.tsv ${prefix}_unaligned.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args      = task.ext.args   ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    def features  = gff             ? "--features $gff" : ''
    def reference = fasta           ? "-r $fasta" : ''

    """
    mkdir -p $prefix/combined_reference
    touch $prefix/report.html
    touch $prefix/metaquast.log
    touch $prefix/combined_reference/report.tsv
    touch $prefix/combined_reference/report.pdf
    touch $prefix/combined_reference/transposed_report.txt
    touch $prefix/combined_reference/transposed_report.tex
    touch $prefix/combined_reference/icarus.html
    touch $prefix/combined_reference/report.tex
    touch $prefix/combined_reference/report.txt

    mkdir -p $prefix/basic_stats
    touch $prefix/combined_reference/basic_stats/cumulative_plot.pdf
    touch $prefix/combined_reference/basic_stats/Nx_plot.pdf
    touch $prefix/combined_reference/basic_stats/genome_GC_content_plot.pdf
    touch $prefix/combined_reference/basic_stats/GC_content_plot.pdf

    mkdir -p $prefix/icarus_viewers
    touch $prefix/combined_reference/icarus_viewers/contig_size_viewer.html

    ln -s $prefix/combined_reference/report.tsv ${prefix}.tsv

    if [ $fasta ]; then
        touch $prefix/combined_reference/basic_stats/NGx_plot.pdf
        touch $prefix/combined_reference/basic_stats/gc.icarus.txt

        mkdir -p $prefix/combined_reference/aligned_stats
        touch $prefix/combined_reference/aligned_stats/NAx_plot.pdf
        touch $prefix/combined_reference/aligned_stats/NGAx_plot.pdf
        touch $prefix/combined_reference/aligned_stats/cumulative_plot.pdf

        mkdir -p $prefix/combined_reference/contigs_reports
        touch $prefix/combined_reference/contigs_reports/all_alignments_transcriptome.tsv
        touch $prefix/combined_reference/contigs_reports/contigs_report_transcriptome.mis_contigs.info
        touch $prefix/combined_reference/contigs_reports/contigs_report_transcriptome.stderr
        touch $prefix/combined_reference/contigs_reports/contigs_report_transcriptome.stdout
        touch $prefix/combined_reference/contigs_reports/contigs_report_transcriptome.unaligned.info
        mkdir -p $prefix/combined_reference/contigs_reports/minimap_output
        touch $prefix/combined_reference/contigs_reports/minimap_output/transcriptome.coords
        touch $prefix/combined_reference/contigs_reports/minimap_output/transcriptome.coords.filtered
        touch $prefix/combined_reference/contigs_reports/minimap_output/transcriptome.coords_tmp
        touch $prefix/combined_reference/contigs_reports/minimap_output/transcriptome.sf
        touch $prefix/combined_reference/contigs_reports/minimap_output/transcriptome.unaligned
        touch $prefix/combined_reference/contigs_reports/minimap_output/transcriptome.used_snps
        touch $prefix/combined_reference/contigs_reports/misassemblies_frcurve_plot.pdf
        touch $prefix/combined_reference/contigs_reports/misassemblies_plot.pdf
        touch $prefix/combined_reference/contigs_reports/misassemblies_report.tex
        touch $prefix/combined_reference/contigs_reports/misassemblies_report.tsv
        touch $prefix/combined_reference/contigs_reports/misassemblies_report.txt
        touch $prefix/combined_reference/contigs_reports/transcriptome.mis_contigs.fa
        touch $prefix/combined_reference/contigs_reports/transposed_report_misassemblies.tex
        touch $prefix/combined_reference/contigs_reports/transposed_report_misassemblies.tsv
        touch $prefix/combined_reference/contigs_reports/transposed_report_misassemblies.txt
        touch $prefix/combined_reference/contigs_reports/unaligned_report.tex
        touch $prefix/combined_reference/contigs_reports/unaligned_report.tsv
        touch $prefix/combined_reference/contigs_reports/unaligned_report.txt

        mkdir -p $prefix/combined_reference/genome_stats
        touch $prefix/combined_reference/genome_stats/genome_info.txt
        touch $prefix/combined_reference/genome_stats/transcriptome_gaps.txt
        touch $prefix/combined_reference/icarus_viewers/alignment_viewer.html

        ln -sf ${prefix}/combined_reference/contigs_reports/misassemblies_report.tsv ${prefix}_misassemblies.tsv
        ln -sf ${prefix}/combined_reference/contigs_reports/unaligned_report.tsv ${prefix}_unaligned.tsv
        ln -sf ${prefix}/combined_reference/contigs_reports/all_alignments_transcriptome.tsv ${prefix}_transcriptome.tsv

    fi

    if ([ $fasta ] && [ $gff ]); then
        touch $prefix/combined_reference/genome_stats/features_cumulative_plot.pdf
        touch $prefix/combined_reference/genome_stats/features_frcurve_plot.pdf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
