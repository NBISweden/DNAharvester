#! /usr/bin/env nextflow

process SAMREMOVEDUP {
    tag "$meta.id"
    label 'process_samremovedup'

    conda "conda-forge::python=3.13.0 bioconda::samtools=1.21"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools_python:3740ccb1c28d345d' :
        'community.wave.seqera.io/library/samtools_python:44d9e6d118af41a9' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*dedup.bam")         , emit: dedup
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref_prefix = task.ext.ref_prefix ?: "${meta2.id}".replaceAll(/\.(fasta|fna|fa)$/, '')

    """
    samtools view \\
        -@ ${task.cpus-1} \\
        -h $bam | \\
        samremovedup.py | \\
        samtools view \\
        -b \\
        -o ${prefix}.${ref_prefix}.dedup.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
