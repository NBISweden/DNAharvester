process SAMTOOLS_FAIDX {
    tag "$fasta"
    label 'process_single'

    conda "bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools:1.20--ad906e74fde1812b' :
        'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464' }"

    input:
    path(fasta)

    output:
    path ("*.fai")                         , emit: fai, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    if [[ $fasta == *.gz ]]; then
        gunzip -c $fasta > ${fasta%.gz}
    fi

    samtools \\
        faidx \\
        \${fasta%.gz} \\
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
    touch \${fasta%.gz}.fai

    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
