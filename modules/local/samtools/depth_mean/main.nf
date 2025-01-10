process SAMTOOLS_DEPTH_MEAN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::htslib=1.21 bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9edc2564215d5cd137a8b25ca8a311600987186d406b092022444adf3c4447f7/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.21--6cb89bfd40cbaabf' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(intervals)

    output:
    tuple val(meta), path("*.dpstats.txt"), emit: dpstats
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def positions = intervals ? "-b ${intervals}" : ""
    """
    samtools \\
        depth \\
        --threads ${task.cpus-1} \\
        $args \\
        $positions \\
        -o ${prefix}.tsv \\
        $bam

    awk \\
        '{sum+=\$3} END { print sum/NR }' \\
        ${prefix}.tsv | \\
        awk \\
        '{ printf "%.0f", \$1 }' \\
        > ${prefix}.dpstats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        awk: \$(echo \$(awk --version 2>&1) | sed 's/^.*awk //; s/Using.*\$//')
    END_VERSIONS
    """
}
