process PLOT_PRESEQ {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::matplotlib=3.8.4"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib:df831ede8509d5c2' :
        'community.wave.seqera.io/library/matplotlib:3.8.4--3dafe3963199ad72' }"

    input:
    tuple val(meta), path(preseq_txt)

    output:
    tuple val(meta), path("*_preseq.png")       , emit: preseq_plot
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in {{ name }}/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plot_preseq.py \\
        ${preseq_txt} \\
        ${args} \\
        ${prefix}_preseq.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}