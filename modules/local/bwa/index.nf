process BWA_INDEX {
    tag "$fasta"
    label 'process_medium'

    conda "bioconda::bwa=0.7.18 bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bwa_samtools:813d9b5fea3890ec' :
        'community.wave.seqera.io/library/bwa_samtools:3938c84206f62975' }"

    storeDir "${output_dir}"

    input:
    tuple val(meta2), path(fasta)
    val (output_dir)

    output:
    tuple val(meta2), path("${meta2.id}.{amb,ann,bwt,pac,sa}", arity: '5')  , emit: index
    tuple val(meta2), val(output_dir)                                       , emit: index_dir
    path "versions.yml"                                                     , emit: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bwa \\
        index \\
        $args \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}