#! /usr/bin/env nextflow

process SEQ_STATS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(json)
    tuple val(meta), path(raw_bam_flagstat)
    tuple val(meta), path(dedup_lib_flagstat)

    output:
    tuple val(meta), path("*uniq_id.stats.txt"), emit: stats_txt

    script:

    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ ! -f uniq_id.stats.txt ]; then
        printf "sample_ID\\ttotal reads\\tmapped_reads\\tuniq_reads\\n" > uniq_id.stats.txt
    fi
    printf "${prefix}\\t" >> uniq_id.stats.txt
    cat $json | grep -A 1 "read1_before_filtering" | awk -F ': ' '/"total_reads"/ {gsub(/,/, "", \$2); printf \$2 "\\t"}' >> uniq_id.stats.txt
    cat $raw_bam_flagstat | grep -m 1 "mapped (" | awk '{printf \$1 "\\t"}' >> uniq_id.stats.txt
    cat $dedup_lib_flagstat | grep -m 1 "mapped (" | awk '{printf \$1 "\\n"}' >> uniq_id.stats.txt










    """
}
