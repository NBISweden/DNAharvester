process BCFTOOLS_VARIANT_CALLING {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h3a4d415_1' :
        'quay.io/biocontainers/bcftools:1.21--h3a4d415_1' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("*_sorted.bcf")              , emit: bcftools_sorted_bcf
    tuple val(meta), path("*.bcf.stats")               , emit: bcftools_bcf_stats
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mapQ = task.ext.mapQ ?: "${params.mapping_quality}"
    def baseQ = task.ext.baseQ ?: "${params.base_quality}"
    def non_variant_sites = params.keep_no_variant_sites ? '' : '-v'


    """
    bcftools mpileup \\
        -q ${mapQ} \\
        -Q ${baseQ} \\
        -B \\
        -Ou \\
        -f ${fasta} \\
        ${bam} \\
        --ignore-RG \\
        $args \\
        --threads ${task.cpus-1} | \\
    bcftools call \\
        -m \\
        ${non_variant_sites} \\
        -Ob \\
        --threads ${task.cpus-1} \\
        -o ${prefix}.bcf


    rm ${prefix}.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version | head -n1 | awk '{print \$2}')"
    END_VERSIONS
    """
}

