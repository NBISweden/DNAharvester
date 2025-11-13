process CREATE_AMBER_SAMPLESHEET {
    tag "$meta.id"
    label 'process_create_amber_samplesheet'

    conda "conda-forge::python=3.13.0"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.13.0--a8086dc1de1c4e39' :
        'community.wave.seqera.io/library/python:3.13.0--a025ad9838d75455' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.amber.tsv") , emit: tsv
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in {{ name }}/bin/
    """
    create_amber_samplesheet.py \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}