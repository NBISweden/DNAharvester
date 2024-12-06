process CONCATENATE_TARGET_DECOY_REFS {
    label 'process_single'

    input:
    tuple path(fasta), path(decoy)

    output:
    path "*.fasta"     , emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat ${fasta} ${decoy} > ${fasta}-${decoy}.fasta
    """
}