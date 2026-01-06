process FILTERBAM_PLOT {
    tag "$meta.id"
    label 'process_filterbam_plot'

    conda "conda-forge::matplotlib=3.8.4"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_matplotlib_numpy_python_pruned:84686eb793124fcd' :
        'community.wave.seqera.io/library/pysam_matplotlib_numpy_python_pruned:7a2de054bdadda21' }"

    input:
    tuple val(meta), path(bam), path(bai), path(filterBAM_table)


    output:
    tuple val(meta), path("*_ms_plot.pdf")       , emit: ms_plot
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in {{ name }}/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    microbial_screening_plot.py \\
        --bam_file ${bam} \\
        --out ${prefix}_ms_plot.pdf \\
        --filterBAM_table ${filterBAM_table}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}



