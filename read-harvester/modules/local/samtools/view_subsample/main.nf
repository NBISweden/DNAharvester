process SAMTOOLS_VIEW_SUBSAMPLE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools:1.20--ad906e74fde1812b' :
        'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464' }"

    input:
    tuple val(meta), path(bam)
    path(fasta)

    output:
    tuple val(meta), path("*.bam") , emit: subsampled_bam
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def subsample = params.subsample ?: '1'
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    total=\$(samtools view -F 4 -q 1 -c $bam) # total number of reads, not pairs (excluding unmapped and duplicated multi-aligned reads)
    frac=\$(awk -v s=$subsample -v t=\$total "BEGIN {print s/t}") # fraction of reads to keep

    samtools \\
        view \\
        --threads ${task.cpus-1} \\
        ${reference} \\
        -s \$frac \\
        $args \\
        -o ${prefix}.approx_${subsample}_reads.bam \\
        $bam \\
        $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def subsample = task.ext.subsample ?: '1'
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.approx_${subsample}_reads.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}