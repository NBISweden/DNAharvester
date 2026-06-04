process BCFTOOLS_CALL {
    tag "$meta.id"
    label 'process_bcftools_call'

    conda "bioconda::bcftools=1.21"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h3a4d415_1' :
        'quay.io/biocontainers/bcftools:1.21--h3a4d415_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta_fasta), path(fasta)
    tuple val(meta_fai), path(fai)
    tuple val(meta2), path(bed_file)

    output:
    tuple val(meta), path("*_sorted.bcf")              , emit: sorted_bcf
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mapQ = task.ext.mapQ ?: "${params.mapping_quality}"
    def baseQ = task.ext.baseQ ?: "${params.bcftools_base_quality}"
    def non_variant_sites = task.ext.non_variant_sites ?: params.bcftools_keep_non_variant_sites.toBoolean() ? '-v' : ''
    def regions_file = (meta2.id == 'null' && bed_file.name == 'null') ? "" : "-R ${bed_file}"


    """
    bcftools mpileup \\
        -q ${mapQ} \\
        -Q ${baseQ} \\
        -Ou \\
        -f ${fasta} \\
        -a FORMAT/DP \\
        ${bam} \\
        --ignore-RG \\
        ${regions_file} \\
        $args \\
        --threads ${task.cpus-1} | \\
    bcftools call \\
        -m \\
        ${non_variant_sites} \\
        -Ob \\
        --threads ${task.cpus-1} \\
        -o ${prefix}.bcf

    bcftools sort -O b -o ${prefix}_sorted.bcf ${prefix}.bcf
    bcftools index ${prefix}_sorted.bcf

    rm ${prefix}.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version | head -n1 | awk '{print \$2}')"
    END_VERSIONS
    """
}

