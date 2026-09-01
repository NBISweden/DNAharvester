process MAPPING_METRICS_SAMPLE {
    tag "$meta.id"
    label 'process_mapping_metrics_sample'

    conda "bioconda::samtools=1.21 conda-forge::gawk=5.3.1 conda-forge::jq=1.8.1"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools_gawk_jq:0c240a96dd50ca44' :
        'community.wave.seqera.io/library/samtools_gawk_jq:d2a006aa774b3346' }"



    input:
    tuple val(meta),
    path(fastp_json),
    path(raw_bam_flagstat),
    path(filtered_bam_flagstat),
    path(dedup_lib_flagstat),
    path(dedup_lib),
    path(dpstats),
    path(decoy_flagstat)


    output:
    tuple val(meta), path("*.stats.tsv")    , emit: stats_tsv
    path "versions.yml"                     , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = task.ext.prefix ?: params.reference.split('/')[-1].replaceFirst(/\.fasta$|\.fa$/, '')
    def mapping_tool = task.ext.mapping_tool ?: (meta.sample_type == 'ancient' ? params.mapping_tool_ancient : params.mapping_tool_modern)

    """

    ## collect stats
    id="${prefix}"
    raw_reads_pairs=\$(cat ${fastp_json} | jq '.read1_before_filtering.total_reads' | awk '{sum += \$1} END {print sum}')
    processed_reads=\$(awk '\$NF=="primary"{sum+=\$1} END{print sum+0}' ${raw_bam_flagstat})
    reference=\$(basename ${reference})
    mapping_tool="${mapping_tool}"
    mapped_reads=\$(cat ${raw_bam_flagstat} | grep "primary mapped (" | awk '{sum += \$1} END {print sum}')
    decoy_reads=\$(cat ${decoy_flagstat} | grep "primary mapped (" | awk '{sum += \$1} END {print sum}')
    mq_filter=${params.mapping_quality}
    filtered_reads=\$(cat ${filtered_bam_flagstat} | grep "primary mapped (" | awk '{sum += \$1} END {print sum}')
    filtered_paired_reads=\$(awk '/paired in sequencing/{sum+=\$1} END{print sum+0}' ${filtered_bam_flagstat})
    mapped_fragments=\$(awk -v filtered_reads=\$filtered_reads -v filtered_paired_reads=\$filtered_paired_reads 'BEGIN { print filtered_reads - (filtered_paired_reads/2) }')
    endogenous_DNA=\$(awk -v mapped_fragments=\$mapped_fragments -v raw_reads_pairs=\$raw_reads_pairs 'BEGIN { if (raw_reads_pairs > 0) print ((mapped_fragments / raw_reads_pairs)*100); else print 0 }')
    uniq_reads=\$(cat ${dedup_lib_flagstat} | grep "primary mapped (" | awk '{sum += \$1} END {print sum}')
    library_complexity=\$(awk -v uniq_reads=\$uniq_reads -v filtered_reads=\$filtered_reads 'BEGIN { if (filtered_reads > 0) print ((uniq_reads / filtered_reads)*100); else print 0 }')
    samtools stats --threads ${task.cpus} ${dedup_lib} > ${prefix}-samtools-stats
    depth=\$(cat ${dpstats})
    min_reads_len=\$(awk '/^RL/ {print \$2}' ${prefix}-samtools-stats | head -n 1)
    max_reads_len=\$(awk '/^SN/ && /maximum length/ {print \$4}' ${prefix}-samtools-stats)
    mean_reads_len=\$(awk '/^SN/ && /average length/ {print \$4}' ${prefix}-samtools-stats)
    median_reads_len=\$(awk '/^RL/ {total+=\$3; lengths[\$2]=\$3} END {median=total/2; sum=0; for (len in lengths) {sum+=lengths[len]; if (sum>=median) {print len; break}}}' ${prefix}-samtools-stats)

    ## write header and row
    HEADER="id\\traw_reads_pairs\\tprocessed_reads\\tref_genome\\tmapping_tool\\tmapped_reads"
    ROW="\$id\\t\$raw_reads_pairs\\t\$processed_reads\\t\$reference\\t\$mapping_tool\\t\$mapped_reads"

    ## add optional decoy

    if [[ "${decoy_flagstat}" != "null" && "${decoy_flagstat}" != "/dev/null" ]]; then
        HEADER+="\\tdecoy_reads"
        decoy_reads=\$(cat ${decoy_flagstat} | grep "primary mapped (" | awk '{sum += \$1} END {print sum}')
        ROW+="\\t\$decoy_reads"
    fi

    HEADER+="\\tmq_filter\\tfiltered_reads\\tendogenous_DNA\\tunique_reads\\tlibrary_complexity\\tdepth\\tmin_read_length\\tmax_read_length\\tmean_read_length\\tmedian_read_length"
    ROW+="\\t\$mq_filter\\t\$filtered_reads\\t\$endogenous_DNA\\t\$uniq_reads\\t\$library_complexity\\t\$depth\\t\$min_reads_len\\t\$max_reads_len\\t\$mean_reads_len\\t\$median_reads_len"

    ## write output
    printf "\$HEADER\\n" > ${prefix}.stats.tsv
    printf "\$ROW\\n" >> ${prefix}.stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        awk: \$(awk --version | head -n 1 | awk '{print \$1, \$2, \$3}')
        jq: \$(jq --version)
    END_VERSIONS
    """
}
