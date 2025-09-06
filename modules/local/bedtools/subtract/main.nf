process BEDTOOLS_SUBTRACT {
    label 'process_low'

    conda "bioconda::bedtools=2.31.1"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta2), path(intervals1)
    tuple val(meta3), path(intervals2)

    output:
    tuple val(meta2), path("*.bed") , emit: bed
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix1 = task.ext.prefix1 ?: "${meta2.id}"
    def prefix2 = task.ext.prefix2 ?: "${meta3.id}"
    """
    bedtools \\
        subtract \\
        -a $intervals1 \\
        -b $intervals2 \\
        $args \\
        > ${prefix1}.subtract-${prefix2}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}