process KRAKEN2 {
    tag "$meta.id"
    label 'process_kraken2'

    conda "bioconda::kraken2=2.1.3"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/kraken2:eed6d8ea184673ff' :
        'community.wave.seqera.io/library/kraken2:3773f4955380979e' }"

    input:
    tuple val(meta) , path(reads)
    path(kraken2_database)

    output:
    tuple val(meta), path("*.kraken2")     , emit: kraken2_output
    tuple val(meta), path("*.output")      , emit: kraken2_report
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    kraken2 \\
        $reads $args \\
        --db ${kraken2_database} \\
        --threads ${task.cpus} \\
        --report-minimizer-data \\
        --use-names \\
        --output ${prefix}.kraken2 \\
        --report ${prefix}.output


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Version: //; s/ .*//')
    END_VERSIONS
    """
}