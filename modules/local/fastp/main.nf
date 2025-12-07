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
    // This channel will emit:
    // - [meta, [read.fastq.gz]] for single-end
    // - [meta, [merged.fastq.gz]] for paired-end (merged)
    // - [meta, [R1.fastq.gz, R2.fastq.gz]] for paired-end (unmerged)
    tuple val(meta), path('*.fastp.fastq.gz')           , emit: processed_reads
    tuple val(meta), path('*.json')                     , emit: json
    tuple val(meta), path('*.html')                     , emit: html
    tuple val(meta), path('*.log')                      , emit: log
    // Optional outputs for unmerged reads when merging is performed
    tuple val(meta), path('*.unmerged.R*.fastq.gz')     , optional:true, emit: unmerged_reads
    path "versions.yml"                                 , emit: versions

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
        def merge_cmd = params.merge_reads.toBoolean() ?
            "--merge --merged_out ${prefix}.merged.fastp.fastq.gz --correction --overlap_len_require 15 --overlap_diff_limit 1 --out1 ${prefix}.unmerged.R1.fastq.gz --out2 ${prefix}.unmerged.R2.fastq.gz" :
            "--out1 ${prefix}_R1.fastp.fastq.gz --out2 ${prefix}_R2.fastp.fastq.gz"
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
    fastp \\
        --in1 ${prefix}_1.fastq.gz \\
        --in2 ${prefix}_2.fastq.gz \\
        $merge_cmd \\
        --detect_adapter_for_pe \\
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
