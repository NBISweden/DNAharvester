process BCFTOOLS_RM_INDELS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h3a4d415_1' :
        'quay.io/biocontainers/bcftools:1.21--h3a4d415_1' }"

    input:
    tuple val(meta), path(bcf)

    output:
    tuple val(meta), path("*_rm-indels.bcf") , emit: rm_indels_bcf
    path "versions.yml"                      , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bcftools filter \\
        -i 'INDEL=0' \\
        -Oz \\
        -o ${bcf.baseName}_rm-indels.bcf \\
        ${args} \\
        --threads ${task.cpus-1} \\
        ${bcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version | head -n1 | awk '{print \$2}')"
    END_VERSIONS
    """
}