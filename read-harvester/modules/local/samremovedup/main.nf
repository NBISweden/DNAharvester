#! /usr/bin/env nextflow

process SAMREMOVEDUP {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.11.0 bioconda::samtools=1.16.1"
    container "quay.io/biocontainers/mulled-v2-1a35167f7a491c7086c13835aaa74b39f1f43979:9254eac8981f615fb6c417fa44e77c3b44bc3abd-0"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*dedup.bam")                       , emit: dedup
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view \\
        -@ ${task.cpus-1} \\
        -h $bam | \\
        python3 samremovedup.py | \\
        samtools view \\
        -b \\
        -o ${prefix}.dedup.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
