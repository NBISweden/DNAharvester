process SAMTOOLS_DEPTH {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.20 bioconda::htslib=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/htslib_samtools:1.20--9fb9031594b6902c' :
        'community.wave.seqera.io/library/htslib_samtools:1.20--11a4e6daa46930ec' }"

    input:
    tuple val(meta), path(bam)
    path(intervals)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
