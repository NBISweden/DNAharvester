process FAI_TO_BED {
    tag "$meta2.id"
    label 'process_single'

    input:
    tuple val(meta2), path(fai)

    output:
    tuple val(meta2), path "*.bed"     , emit: bed

    when:
    task.ext.when == null || task.ext.when
    def prefix = task.ext.prefix ?: "${meta2.id}"

    script:
    """
    awk -v OFS='\t' '{print \$1, "0", \$2}' ${fai} > ${prefix}.bed
    """
}