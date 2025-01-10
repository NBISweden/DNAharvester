#! /usr/bin/env nextflow

process SAMREMOVEDUP {
    tag "$meta.id"
    label 'process_high'

    conda "conda-forge::python=3.13.0 bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools_python:3740ccb1c28d345d' :
        'community.wave.seqera.io/library/samtools_python:44d9e6d118af41a9' }"

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
