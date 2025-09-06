process CONSENSUS_CALL_MIA {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.13.0"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.13.0--a8086dc1de1c4e39' :
        'community.wave.seqera.io/library/python:3.13.0--a025ad9838d75455' }"

    input:
    tuple val(meta), path(mia_maln_41)

    output:
    tuple val(meta), path("*.fasta")    , emit: mia_fasta
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in {{ name }}/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    def coverage_threshold = task.ext.coverage_threshold ?: "${params.coverage_threshold}"
    def call_fraction = task.ext.call_fraction ?: "${params.call_fraction}"
    def quality_threshold = task.ext.quality_threshold ?: "${params.quality_threshold}"

    """
    consensus_call_mia.py \\
        -c ${coverage_threshold} \\
        -p ${call_fraction} \\
        -q ${quality_threshold} \\
        -I ${prefix} \\
        -m ${mia_maln_41} \\
        -o ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}