process RM_SHORT_READS {
    tag "$meta.id"
    label 'process_rm_short_reads'

    conda "bioconda::htslib=1.21 bioconda::samtools=1.21"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9edc2564215d5cd137a8b25ca8a311600987186d406b092022444adf3c4447f7/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.21--6cb89bfd40cbaabf' }"

    input:
    tuple val(meta), path(bam), path(read_len)

    output:
    tuple val(meta), path("*.bam")          , emit: bam
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    len=\$(head -n 1 $read_len | awk -F ': ' '{print \$2}')
    samtools view -h --threads ${task.cpus} $bam | \\
    awk -v len=\$len 'length(\$10) >= len || \$1 ~ /^@/' | \\
    samtools view -bh --threads ${task.cpus} -o ${prefix}-rl\${len}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        awk: \$(echo \$(awk --version 2>&1) | sed 's/^.*awk //; s/Using.*\$//')
    END_VERSIONS
    """
}