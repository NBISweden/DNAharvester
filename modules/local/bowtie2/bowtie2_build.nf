process BOWTIE2_BUILD {
    tag "$fasta"
    label 'process_bowtie2_build'

    conda "bioconda::bowtie2=2.5.4"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bowtie2:2.5.4--2ec535d45cd82f0b' :
        'quay.io/biocontainers/bowtie2:2.5.4--he96a11b_6' }"

    storeDir "${output_dir}"

    input:
    tuple val(meta2), path(fasta)
    val (output_dir)

    output:
    tuple val(meta2), path("${meta2.id}*bt2", arity: '6')       , emit: index
    tuple val(meta2), val(output_dir)                           , emit: index_dir
    path "versions.yml"                                         , emit: versions, optional: true
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bowtie2-build \\
        ${args} \\
        --threads ${task.cpus} \\
        ${fasta} ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | grep -i version | head -n 1 | sed 's/.*version //')
    END_VERSIONS
    """
}