process BCFTOOLS_RM_ALLELIC_IMBALANCE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h3a4d415_1' :
        'quay.io/biocontainers/bcftools:1.21--h3a4d415_1' }"

    input:
    tuple val(meta), path(bcf)

    output:
    tuple val(meta), path("*_rmindels.bcf") , emit: rmindels_bcf
    path "versions.yml"                     , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:

    """
    bcftools view --threads ${task.cpus}-e 'GT="0/1" & (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]) < 0.2' ${bcf} | \\
    bcftools view --threads ${threads} -e 'GT="0/1" & (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]) > 0.8' | \\
    bcftools view --threads ${threads} -e 'GT="1/2" & (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]) < 0.2' | \\
    bcftools view --threads ${threads} -e 'GT="1/2" & (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]) > 0.8' \\
    -Ob -o ${bcf.baseName}_rmindels.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version | head -n1 | awk '{print \$2}')"
    END_VERSIONS
    """
}