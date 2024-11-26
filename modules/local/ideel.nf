process IDEEL {
    label 'process_single'

    container "quay.io/microbiome-informatics/ideel:1.0.0"

    input:
    tuple val(meta), val(stats)
    val(before_after)
    
    output:
    tuple val(meta), env(platform),     emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out = "${prefix}_${before_after}.png"
    """
    ideel.py ${stats} ${out}

    """
}