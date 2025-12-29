process SORT_METRICS {
    tag "${input_file.getBaseName()}"
    label 'process_sort_metrics'

    input:
    path(input_file)

    output:
    path("*sorted.tsv"), emit: sorted_file

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: input_file.getBaseName()

    """
    head -n 1 ${input_file} > ${prefix}.sorted.tsv
    tail -n +2 ${input_file} | sort -k1,1  >> ${prefix}.sorted.tsv

    """
}
