process SORT_IDXSTATS {
    tag "$meta.id"
    label 'process_sort_idxstats'

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("*.sorted.tsv"),  emit: sorted_idxstats

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sort -t'\t' -k3,3nr ${input_file} > ${prefix}.idxstats.sorted.tsv

    """
}