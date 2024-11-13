process HAPLOTOFASTA {

    conda "conda-forge::python=3.13.0 conda-forge::pandas=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pandas_python:b56517cc205d1f0a' :
        'community.wave.seqera.io/library/pandas_python:fd8290c2da2fd6ae' }"

    input:
    tuple val(meta), path(haplo)
    path(fai)

    output:
    tuple val(meta), path("*.haplo.fasta") , emit: fasta
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in {{ name }}/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    haplo2fasta.py \\
        ${prefix}.haplo.gz \\
        ${fai} \\
        > ${prefix}.haplo.fasta &&

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}