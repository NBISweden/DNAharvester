process BCFTOOLS_FILTER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h3a4d415_1' :
        'quay.io/biocontainers/bcftools:1.21--h3a4d415_1' }"

    input:
    tuple val(meta), path(sorted_bcf)

    output:
    tuple val(meta), path("*.bcf")      , emit: filtered_bcf
    path "versions.yml"                 , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def variant_quality = task.ext.variant_quality ?: "${params.bcftools_variant_quality}"
    def variant_min_depth = task.ext.variant_min_depth ?: "${params.bcftools_variant_min_depth}"
    def variant_max_depth = task.ext.variant_max_depth ?: "${params.bcftools_variant_max_depth}"
    def variant_gap_indels = task.ext.variant_gap_indels ?: "${params.bcftools_variant_gap_indels}"


    """
    bcftools filter \\
        -e 'QUAL<${variant_quality} || FMT/DP<${variant_min_depth} || FMT/DP>${variant_max_depth}' \\
        -g ${variant_gap_indels} \\
        -Oz \\
        -o ${prefix}_sorted_qual${variant_quality}_dp${variant_min_depth}-${variant_max_depth}_gapindels${variant_gap_indels}.bcf \\
        ${args} \\
        --threads ${task.cpus-1} \\
        ${sorted_bcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version | head -n1 | awk '{print \$2}')"
    END_VERSIONS
    """
}
