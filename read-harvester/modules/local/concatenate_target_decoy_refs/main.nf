process CONCATENATE_TARGET_DECOY_REFS {
    label 'process_single'

    input:
    tuple path(fasta), path(decoy)

    output:
    path "*.fasta"     , emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if [[ ${fasta} != *.gz ]]; then
        gzip -c ${fasta} > ${fasta}.gz
        fasta=${fasta}.gz
    fi

    if [[ ${decoy} != *.gz ]]; then
        gzip -c ${decoy} > ${decoy}.gz
        decoy=${decoy}.gz
    fi

    zcat ${fasta} ${decoy} > concatenated.fasta
    """
}