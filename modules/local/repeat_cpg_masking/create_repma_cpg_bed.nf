process CREATE_REPMA_CPG_BED {
    tag "$fasta"
    label 'process_low'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(genomefile)
    tuple val(meta2), path(sorted_repeats_bed)
    tuple val(meta2), path(ref_bed)
    tuple val(meta2), path(cpg_bed)

    output:
    tuple val(meta2), path("*noCpG_ref.bed")                , emit: no_cpg_bed
    tuple val(meta2), path("*CpG_ref.repeats.sorted.bed")   , emit: cpg_repeats_sorted_bed
    tuple val(meta2), path("*noCpG_ref.repma.bed")          , emit: no_cpg_repma_bed
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta2.id}"
    def args = task.ext.args ?: ''

    """
    bedtools subtract -a ${ref_bed} -b ${cpg_bed} > ${prefix}.noCpG_ref.bed
    cat ${sorted_repeats_bed} ${cpg_bed} | sort -k1,1 -k2,2n | bedtools merge > ${prefix}.CpG_ref.repeats.bed
    bedtools sort -g ${genomefile} -i ${prefix}.CpG_ref.repeats.bed > ${prefix}.CpG_ref.repeats.sorted.bed
    bedtools subtract -a ${ref_bed} -b ${prefix}.CpG_ref.repeats.sorted.bed > ${prefix}.noCpG_ref.repma.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools -v 2>&1) | awk '{print \$3}')
    END_VERSIONS
    """
}