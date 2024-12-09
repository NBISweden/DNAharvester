process CONCATENATE_TARGET_DECOY_REFS {
    tag "$fasta"
    label 'process_single'

    input:
    path fasta
    path decoy

    output:
    path "concatenated.fasta", emit: concatenated

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if [[ ${fasta} == *.gz ]]; then
        gunzip -c ${fasta} > ${fasta%.gz}
        fasta=${fasta%.gz}
    fi

    if [[ ${decoy} == *.gz ]]; then
        gunzip -c ${decoy} > ${decoy%.gz}
        decoy=${decoy%.gz}
    fi

    cat ${fasta} ${decoy} > concatenated.fasta
    """
}