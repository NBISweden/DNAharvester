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
    # Split the *.haplo.gz file per chromosome
    zcat < ${prefix}.haplo.gz | \\
        awk 'FNR == 1 {next}{print>(\$1".haplo")}' &&

    # Loop through the chromosome *.haplo files and convert them to fasta format
    # to avoid storing the genome-wide *.haplo.gz file into memory
    for chr in \$(ls *.haplo); do 
        gzip \${chr} &&
        haplo2fasta.py \\
            \${chr}.gz \\
            ${fai} \\
            > \${chr}.fasta
    done &&

    # Concatenate the chromosome *.haplo.fasta files 
    cat *.haplo.fasta > ${prefix}.haplo.fasta &&

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}