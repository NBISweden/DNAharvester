#! /usr/bin/env nextflow

process SAMREMOVEDUP {
    tag "$meta.id"
    label 'process_high'

    conda "conda-forge::python=3.13.0 bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools_python:e886d0f7a340b459' :
        'community.wave.seqera.io/library/samtools_python:97109fdca4337830' }"

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
        samremovedup.py | \\
        samtools view \\
        -b \\
        -o ${prefix}.dedup.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
