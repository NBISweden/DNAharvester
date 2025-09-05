process SAMTOOLS_VIEW_REGIONS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::htslib=1.21 bioconda::samtools=1.21"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9edc2564215d5cd137a8b25ca8a311600987186d406b092022444adf3c4447f7/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.21--6cb89bfd40cbaabf' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(intervals)

    output:
    tuple val(meta), path("*.bam"),                                    emit: bam
    path  "versions.yml",                                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    prefix2 = task.ext.prefix2 ?: "${meta2.id}"
    def positions = intervals ? "-L ${intervals}" : ""
    """
    samtools \\
        view \\
        --threads ${task.cpus-1} \\
        $args \\
        $positions \\
        -o ${prefix}-${prefix2}-regions.bam \\
        $input \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
