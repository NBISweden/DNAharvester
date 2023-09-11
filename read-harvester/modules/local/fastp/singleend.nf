#! /usr/bin/env nextflow

process FASTP_SINGLE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::fastp=0.23.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
        'biocontainers/fastp:0.23.4--h5f740d0_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastp.fastq.gz') , emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    tuple val(meta), path('*.log')            , emit: log
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz
    fastp \\
        --in1 ${prefix}.fastq.gz \\
        --out1 ${prefix}.fastp.fastq.gz \\
        --thread $task.cpus \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        $args \\
        2> ${prefix}.fastp.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}