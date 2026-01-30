process SAMTOOLS_FAIDX {
    tag "$fasta"
    label 'process_samtools_faidx'

    conda "bioconda::htslib=1.21 bioconda::samtools=1.21"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9edc2564215d5cd137a8b25ca8a311600987186d406b092022444adf3c4447f7/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.21--6cb89bfd40cbaabf' }"

    input:
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta2), path ("*.fai")       , emit: fai, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        faidx \\
        $fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def match = (task.ext.args =~ /-o(?:utput)?\s(.*)\s?/).findAll()
    def fastacmd = match[0] ? "touch ${match[0][1]}" : ''
    """
    ${fastacmd}
    touch ${fasta}.fai

    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
