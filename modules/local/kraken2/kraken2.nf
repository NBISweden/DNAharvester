process KRAKEN2 {
    tag "$meta.id"
    label 'process_kraken2'

    conda "bioconda::kraken2=2.1.3 conda-forge::pigz=2.8"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/kraken2_pigz:c88f720548c3d49a' :
        'community.wave.seqera.io/library/kraken2_pigz:be4a80723677f716' }"

    input:
    tuple val(meta) , path(reads)
    path(kraken2_database)

    output:
    tuple val(meta), path("*.kraken2.out")                  , emit: kraken2_output
    tuple val(meta), path("*.report.txt")                   , emit: kraken2_report
    tuple val(meta), path("*.kraken_unclassified.fq.gz")    , emit: kraken2_unclassified
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    kraken2 \\
        --db ${kraken2_database} \\
        --report-minimizer-data \\
        --unclassified-out ${prefix}.kraken_unclassified.fq \\
        --output ${prefix}.kraken2.out \\
        --report ${prefix}.report.txt \\
        --threads ${task.cpus} \\
        ${args} \\
        ${reads}

    pigz -p ${task.cpus} ${prefix}.kraken_unclassified.fq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Version: //; s/ .*//')
    END_VERSIONS
    """
}