process SAMTOOLS_VIEW_SUBSAMPLE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::htslib=1.21 bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9edc2564215d5cd137a8b25ca8a311600987186d406b092022444adf3c4447f7/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.21--6cb89bfd40cbaabf' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.bam") , emit: subsampled_bam
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def subsample = params.subsample ?: '1000000' // default to 1 million reads
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    total=\$(samtools view -F 4 -q 1 -c $bam) # total number of reads, not pairs (excluding unmapped and duplicated multi-aligned reads)
    frac=\$(awk -v s=$subsample -v t=\$total "BEGIN {print s/t}") # fraction of reads to keep

    ## skip subsampling if the total number of reads is less or equal to subsample number
    if (( \$total <= $subsample )); then
        samtools view \\
        --threads \\${task.cpus-1} \\
        -bh $args -F 4 -q 1 \\
        -o ${prefix}.approx_\${total}_reads.bam \\
        $bam
    else
        ## subsample the bam file
        frac=\$(awk -v s=$subsample -v t=\$total "BEGIN {print s/t}") # fraction of reads to keep
        samtools view \\
            --threads ${task.cpus-1} \\
            -s \$frac \\
            -bh $args -F 4 -q 1 \\
            -o ${prefix}.approx_${subsample}_reads.bam \\
            $bam
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}