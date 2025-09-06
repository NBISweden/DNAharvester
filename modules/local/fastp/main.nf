#! /usr/bin/env nextflow
process FASTP {
    tag "$meta.id"
    label 'process_fastp'

    conda "bioconda::fastp=0.24.0"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/fastp:0.24.0--0397de619771c7ae' :
        'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastp.fastq.gz')       , emit: reads //output merged fastq in case of paired-end, output adapter trimmed fastq in case of single-end
    tuple val(meta), path('*.json')                 , emit: json
    tuple val(meta), path('*.html')                 , emit: html
    tuple val(meta), path('*.log')                  , emit: log
    path "versions.yml"                             , emit: versions
    tuple val(meta), path('*.fastp_R1.fastq.gz')    , optional:true, emit: reads_unmerged_R1
    tuple val(meta), path('*.fastp_R2.fastq.gz')    , optional:true, emit: reads_unmerged_R2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
    """
    [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz
    fastp \\
        --in1 ${prefix}.fastq.gz \\
        --out1 ${prefix}.fastp.fastq.gz \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        --thread $task.cpus \\
        $args \\
        2> ${prefix}.fastp.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
    } else {
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
    fastp \\
        --in1 ${prefix}_1.fastq.gz \\
        --in2 ${prefix}_2.fastq.gz \\
        --merge \\
        --merged_out ${prefix}.merged.fastp.fastq.gz \\
        --correction \\
        --overlap_len_require 15 \\
        --overlap_diff_limit 1 \\
        --detect_adapter_for_pe \\
        --out1 ${prefix}.fastp_R1.fastq.gz \\
        --out2 ${prefix}.fastp_R2.fastq.gz \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        --thread $task.cpus \\
        $args \\
        2> ${prefix}.fastp.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
    }
}







