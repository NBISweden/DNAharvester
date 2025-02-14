process FILTER_FASTQ {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::seqtk=1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/seqtk:b3f4cb75b2c12d62' :
        'community.wave.seqera.io/library/seqtk:1.4--f2bbc7882319500b' }"

    input:
    tuple val(meta), path(reads), path(kraken2_output)


    output:
    tuple val(meta), path("*_kraken2-filtered.fastq.gz")    , emit: filtered_reads
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Extract read names where column 1 is "U"
    awk '\$1 == "U" { print \$2 }' $kraken2_output > ${prefix}_unclassified_readnames.txt

    seqtk subseq \\
        $reads \\
        ${prefix}_unclassified_readnames.txt \\
        | gzip > ${prefix}_kraken2-filtered.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*//')
    END_VERSIONS
    """
}