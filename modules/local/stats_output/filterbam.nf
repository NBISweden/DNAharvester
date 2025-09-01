process FILTERBAM {
    tag "$meta.id"
    label 'process_low'

    conda "genomewalker::bam-filter=1.2.1"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cxx-compiler_python_pip_bam-filter:d381468da6310c10' :
        'community.wave.seqera.io/library/cxx-compiler_python_pip_bam-filter:57f86c8e1a5a2597' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_filterBAM.csv"),       emit: filterBAM_stats
    path  "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    filterBAM filter \\
        ${args} \\
        --threads ${task.cpus} \\
        --bam ${bam} \\
        --stats ${bam}_filterBAM.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filterBAM: \$(filterBAM --version 2>&1 | grep filterBAM | head -n 1 | awk '{print \$2}')
    END_VERSIONS

    """
}
