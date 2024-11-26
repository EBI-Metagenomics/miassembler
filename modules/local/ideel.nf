process IDEEL {
    label 'process_single'

    container "quay.io/microbiome-informatics/ideel:1.0.0"

    input:
    tuple val(meta), val(stats)
    
    output:
    tuple val(meta), env(platform),     emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out = "output.png"
    """
    ideel.py ${stats} ${out}


    """
}