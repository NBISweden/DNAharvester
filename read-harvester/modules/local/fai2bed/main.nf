process FAI_TO_BED {
    label 'process_single'

    input:
    path(fai)

    output:
    path "*.bed"     , emit: bed

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk -v OFS='\t' '{print \$1, "0", \$2}' ${fai} > genome.bed
    """
}