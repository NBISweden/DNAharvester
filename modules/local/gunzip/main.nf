process GUNZIP {
    tag "$fasta"
    label 'process_gunzip'

    conda "conda-forge::gzip=1.14"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/gzip:1.14--73d1c92b03f0ea38' :
        'community.wave.seqera.io/library/gzip:1.14--19aaa2c84c85ddbc' }"

    storeDir "${file(params.reference).parent}"

    input:
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta2), path("${fasta.baseName.replaceAll(/\.gz$/, '')}")      , emit: unzip_fasta
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    gzip \\
        -cd \\
        ${args} \\
        ${fasta} \\
        > ${fasta.baseName.replaceAll(/\.gz$/, '')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(echo \$(gzip --version 2>&1) | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """
}