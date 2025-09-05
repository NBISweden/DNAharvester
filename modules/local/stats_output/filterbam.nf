process FILTERBAM {
    tag "$meta.id"
    label 'process_filterbam'

    conda "genomewalker::bam-filter=1.2.1"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cxx-compiler_python_pip_bam-filter:d381468da6310c10' :
        'community.wave.seqera.io/library/cxx-compiler_python_pip_bam-filter:57f86c8e1a5a2597' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*_filterBAM.txt"),       emit: filterBAM_stats
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
        --stats ${prefix}_filterBAM_raw.txt

    ## add sample id to the column 1
    {
        echo -e "ID\t\$(head -n1 ${prefix}_filterBAM_raw.txt)"
        tail -n +2 ${prefix}_filterBAM_raw.txt | sed "s/^/${prefix}\t/"
    } > ${prefix}_filterBAM.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filterBAM: \$(filterBAM --version 2>&1 | grep filterBAM | head -n 1 | awk '{print \$2}')
    END_VERSIONS

    """
}
