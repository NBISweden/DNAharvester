process SEQ_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::samtools=1.21 conda-forge::gawk=5.3.1"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools_gawk:47334f042a11e64b' :
        'community.wave.seqera.io/library/samtools_gawk:2f6ad5ac0a3fef78' }"


    input:
    tuple val(meta),
    path(reads),
    path(fastp_log),
    path(raw_bam_flagstat),
    path(mq_filtered_bam_flagstat),
    path(dedup_lib_flagstat),
    path(dedup_lib)

    output:
    tuple val(meta), path("*.stats.txt")    , emit: stats_txt
    path "versions.yml"                     , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    reference = params.reference.split('/')[-1].replaceFirst(/\.fasta$|\.fa$/, '')

    """
    ## header
    printf "id\\traw_reads\\tmerged_reads\\tref_genome\\tmapping_program\\tmapped_reads\\tmq_filter\\tfiltered_reads\\tunique_reads\\
    \\tmin_read_length\\tmax_read_length\\tmean_read_length\\tmedian_read_lenth\\n" > ${prefix}.stats.txt

    ## collect stats
    id="${prefix}"
    raw_reads=\$(zcat ${reads} | wc -l | awk '{print \$1 / 4}')
    merged_reads=\$(cat ${fastp_log} | grep "Read pairs merged" | awk -F ': ' '{print \$2}')
    reference=\$(basename ${reference})
    mapping_program="bwa aln"
    mapped_reads=\$(cat ${raw_bam_flagstat} | grep -m 1 "mapped (" | awk '{printf \$1 "\\t"}')
    mq_filter=${params.mq}
    filtered_reads=\$(cat ${mq_filtered_bam_flagstat} | grep -m 1 "mapped (" | awk '{printf \$1 "\\t"}')
    uniq_reads=\$(cat ${dedup_lib_flagstat} | grep -m 1 "mapped (" | awk '{printf \$1 "\\n"}')

    samtools stats --threads ${task.cpus} ${dedup_lib} > ${prefix}-samtools-stats
    min_reads_len=\$(awk '/^RL/ {print \$2}' ${prefix}-samtools-stats | head -n 1)
    max_reads_len=\$(awk '/^SN/ && /maximum length/ {print \$4}' ${prefix}-samtools-stats)
    mean_reads_len=\$(awk '/^SN/ && /average length/ {print \$4}' ${prefix}-samtools-stats)
    median_reads_len=\$(awk '/^RL/ {total+=\$3; lengths[\$2]=\$3} END {median=total/2; sum=0; for (len in lengths) {sum+=lengths[len]; if (sum>=median) {print len; break}}}' ${prefix}-samtools-stats)

    ## write stats
    printf "\$id\\t\$raw_reads\\t\$merged_reads\\t\$reference\\t\$mapping_program\\t\$mapped_reads\\t\$mq_filter\\t\$filtered_reads\\t\$uniq_reads\\
    \\t\$min_reads_len\\t\$max_reads_len\\t\$mean_reads_len\\t\$median_reads_len\\n" >> ${prefix}.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        awk: \$(awk --version | head -n 1 | awk '{print \$1, \$2, \$3}')
    END_VERSIONS
    """

}
