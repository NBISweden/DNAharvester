process BWA_INDEX {
    tag "$fasta"
    label 'process_medium'

    conda "bioconda::bwa=0.7.18"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bwa:0.7.18--4543b4091f454101' :
        'community.wave.seqera.io/library/bwa:0.7.18--324359fbc6e00dba' }"

    storeDir "${output_dir}"

    input:
    tuple val(meta2), path(fasta)
    val (output_dir)

    output:
    tuple val(meta2), path("${meta2.id}.{amb,ann,bwt,pac,sa}", arity: '5')  , emit: index
    tuple val(meta2), val(output_dir)                                       , emit: index_dir

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bwa \\
        index \\
        $args \\
        ${fasta}

    """
}