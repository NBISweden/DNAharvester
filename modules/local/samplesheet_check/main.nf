// Modified from https://github.com/nf-core/tools/blob/e5ce6ce20304835bd40f102f038b7e1aadc888b2/nf_core/pipeline-template/modules/local/samplesheet_check.nf

process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_samplesheet_check'

    conda "conda-forge::python=3.13.0"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.13.0--a8086dc1de1c4e39' :
        'community.wave.seqera.io/library/python:3.13.0--a025ad9838d75455' }"

    input:
    path samplesheet

    output:
    path "samplesheet.valid.csv" , emit: csv
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in {{ name }}/bin/
    """
    check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}