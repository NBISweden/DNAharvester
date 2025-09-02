process REFERENCE_DATABASE_UPDATE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::matplotlib=3.8.4"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_matplotlib_numpy_python_pruned:84686eb793124fcd' :
        'community.wave.seqera.io/library/pysam_matplotlib_numpy_python_pruned:7a2de054bdadda21' }"

    input:
    path(accessions_table)
    path(pathogen_reference_database)


    output:
    tuple val(meta), path("*_pathogen_screening_plot.pdf")       , emit: pathogen_screening_plot
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in {{ name }}/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    pathogen_reference_database_update \\
        --table ${accessions_table} \\
        --outpath /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/DNA_HARVESTER/ \\
        --previous_db ${pathogen_reference_database}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}