process RM_SHORT_READS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools:1.20--ad906e74fde1812b' :
        'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam")  , emit: bam
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //#len=$(head -n 1 $read_len | awk -F ': ' '{print \$2}')


    """
    samtools view -h --threads ${task.cpus-1} $bam | \\
    awk 'length(\$10) >= 30 || \$1 ~ /^@/' | \\
    samtools view -bh --threads ${task.cpus-1} -o ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        awk: \$(echo \$(awk --version 2>&1) | sed 's/^.*awk //; s/Using.*\$//')
    END_VERSIONS
    """
}