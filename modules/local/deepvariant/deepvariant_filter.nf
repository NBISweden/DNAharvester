process DEEPVARIANT_FILTER {
    tag "$meta.id"
    label 'process_filter_variants'

    conda "bioconda::bcftools=1.21"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h3a4d415_1' :
        'quay.io/biocontainers/bcftools:1.21--h3a4d415_1' }"

    input:
    tuple val(meta), path(raw_vcf)

    output:
    tuple val(meta), path("*_homalt.bcf")     , path("*_homalt.bcf.csi"), emit: filtered_bcf
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args ?: ''
    def prefix             = task.ext.prefix ?: "${meta.id}"
    def variant_quality    = task.ext.variant_quality    ?: "${params.deepvariant_variant_quality}"
    def variant_min_depth  = task.ext.variant_min_depth  ?: "${params.deepvariant_min_depth}"
    def variant_max_depth  = task.ext.variant_max_depth  ?: "${params.deepvariant_max_depth}"
    def snv_only           = task.ext.snv_only != null ? task.ext.snv_only : params.deepvariant_snv_only.toBoolean()
    def dp                 = "DP${variant_min_depth}-${variant_max_depth}"
    def indel_tag          = snv_only ? '_no_indels' : ''
    def homalt_input       = snv_only ? "${prefix}_${dp}_no_indels.bcf" : "${prefix}_${dp}.bcf"

    """
    ## Filter: depth + quality
    bcftools filter \\
        -i "FMT/DP>=${variant_min_depth} & FMT/DP<${variant_max_depth} & QUAL>=${variant_quality}" \\
        -Ob \\
        --threads ${task.cpus-1} \\
        -o ${prefix}_${dp}.bcf \\
        ${args} \\
        ${raw_vcf}

    ## Conditionally remove indels (SNV-only)
    ${snv_only ? "bcftools view -v snps -Ob --threads ${task.cpus-1} -o ${prefix}_${dp}_no_indels.bcf ${prefix}_${dp}.bcf" : ''}

    ## Keep homozygous-alt only
    bcftools filter \\
        -i 'GT="1/1"' \\
        -Ob \\
        --threads ${task.cpus-1} \\
        -o ${prefix}_${dp}${indel_tag}_homalt.bcf \\
        ${homalt_input}

    bcftools index ${prefix}_${dp}${indel_tag}_homalt.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version | head -n1 | awk '{print \$2}')"
    END_VERSIONS
    """
}
