process PLINK_RECODE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink=1.90b6.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1':
        'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(tfam)
    tuple val(meta), path(tped)

    output:
    tuple val(meta), path("*.vcf.gz")                 , emit: vcf
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plink \\
        --tped ${tped}  \\
        --tfam ${tfam}  \\
        --allow-extra-chr \\
        --threads $task.cpus \\
        --recode vcf-iid bgz \\
        $args \\
        --out $prefix
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(plink --version) | sed 's/^PLINK v//;s/64.*//')
    END_VERSIONS
    """
}