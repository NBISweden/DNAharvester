process CREATE_REPMA_BED {
    tag "$fasta"
    label 'process_low'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(fai)
    tuple val(meta2), path(repeats_bed)

    output:
    tuple val(meta2), path("*.genome")              , emit: genomefile
    tuple val(meta2), path("*_sorted_repeats.bed")  , emit: sorted_repeats_bed
    tuple val(meta2), path("*_ref.bed")             , emit: ref_bed
    tuple val(meta2), path("*_repma.bed")           , emit: repma_bed
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta2.id}"
    def args = task.ext.args ?: ''

    """
    awk -v OFS='\t' '{{print \$1, \$2}}' ${fai} > ${prefix}.genome
    bedtools sort -g ${prefix}.genome -i ${repeats_bed} > ${prefix}_sorted_repeats.bed

    awk -v OFS='\t' '{{print \$1, "0", \$2}}' ${fai} > ${prefix}_ref.bed
    bedtools subtract -a ${prefix}_ref.bed -b ${prefix}_sorted_repeats.bed > ${prefix}_repma.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools -v 2>&1) | awk '{print \$3}')
    END_VERSIONS
    """
}