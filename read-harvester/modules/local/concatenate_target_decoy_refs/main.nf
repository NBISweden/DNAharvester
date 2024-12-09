process CONCATENATE_TARGET_DECOY_REFS {
    tag "$fasta"
    label 'process_single'

    input:
    tuple path(fasta), path(decoy)

    output:
    path "concatenated.fasta", emit: concatenated

    when:
    task.ext.when == null || task.ext.when

    shell:
    """
    if [[ !{fasta} == *.gz ]]; then
        gunzip -c !{fasta} > ${fasta%.gz} &&
        fasta=${fasta%.gz}
    else
        fasta=!{fasta}
    fi

    if [[ !{decoy} == *.gz ]]; then
        gunzip -c !{decoy} > ${decoy%.gz} &&
        decoy=${decoy%.gz}
    else
        decoy=!{decoy}
    fi

    cat ${fasta} ${decoy} > concatenated.fasta
    """
}