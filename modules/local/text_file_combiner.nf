process TXT_COMBINER {
    input:
    tuple val(meta), val(txt1)
    tuple val(meta), val(txt2)

    output:
    tuple val(meta), path("txt_final.txt"), emit: txt_final

    script:
    """
    cat ${txt1} ${txt2} > txt_final.txt
    """
}
