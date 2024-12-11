process ESTIMATE_READ_LEN_CUTOFF {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.13.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.13.0--a8086dc1de1c4e39' :
        'community.wave.seqera.io/library/python:3.13.0--a025ad9838d75455' }"

    input:
    tuple val(meta), path(amber_txt)

    output:
    tuple val(meta), path("*.txt")       , emit: read_len_cutoff
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in {{ name }}/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Define read lenght cutoff threshold 0.05 (5%) for paired-end and 0.01 (1%) for single-end
    def cutoff_threshold = meta.single_end ? 0.01 : 0.05

    """
    select_read_len_cutoff_amber.py \\
        $amber_txt ${cutoff_threshold} > ${prefix}_read_len_cutoff_${cutoff_threshold}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}