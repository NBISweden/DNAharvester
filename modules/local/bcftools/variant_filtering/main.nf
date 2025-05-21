process BCFTOOLS_VARIANT_FILTERING {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h3a4d415_1' :
        'quay.io/biocontainers/bcftools:1.21--h3a4d415_1' }"

    input:
    tuple val(meta), path(bcf)

    output:
    tuple val(meta), path("*_filtered.bcf")           , emit: bcftools_filtered_bcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def qual = task.ext.qual ?: "${params.quality}"


    """
    bcftools filter \\
        -e 'QUAL < ${qual}' \\
        -Oz \\
        -o ${prefix}.filtered.bcf \\
        ${args} \\
        --threads ${task.cpus-1}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version | head -n1 | awk '{print \$2}')"
    END_VERSIONS
    """
}
