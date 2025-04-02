process CREATE_CPG_BED {
    tag "$fasta"
    label 'process_low'

    conda "conda-forge::biopython=1.79"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

    input:
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta2), path("*_cpg.bed")     , emit: cpg_bed
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta2.id}"
    def args = task.ext.args ?: ''

    """
    find_CpG_sites_ref.py \\
        $fasta \\
        ${prefix}_cpg.bed


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | awk '{print \$2}')
    END_VERSIONS
    """
}
