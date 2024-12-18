process PLINK_RECODE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(tfam)
    tuple val(meta), path(tped)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.vcf.gz")                 , emit: vcf
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plink2 \\
        --tfile $prefix  \\
        --allow-extra-chr \\
        --threads $task.cpus \\
        --export vcf id-paste=iid bgz \\
        --ref-from-fa ${fasta} \\
        $args \\
        --out $prefix
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}