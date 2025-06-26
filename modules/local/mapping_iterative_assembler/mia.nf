process MAPPING_ITERATIVE_ASSEMBLER {
    tag "$meta.id"
    label 'process_mia'

    conda "bioconda::mapping-iterative-assembler=1.0"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/mapping-iterative-assembler:1.0--1389c10b012e4570' :
        'quay.io/biocontainers/mapping-iterative-assembler:1.0--h503566f_6' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta4), path(mt_reference)

    output:
    tuple val(meta), path("*.maln.*.41")   , emit: mia_maln_41
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '' // ancient DNA parameters are added via args in the configs/modules.config file
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref_prefix = task.ext.ref_prefix ?: "${meta4.id.replaceAll(/\.(fasta|fna|fa)$/, '')}"

    """
    gunzip -c ${reads} > ${prefix}.unzipped.fastq

    mia -c -C -U -i -F -k 14  \\
        $args \\
        -r ${mt_reference} \\
        -f ${prefix}.unzipped.fastq \\
        -m ${prefix}.${ref_prefix}.maln \\

    ## Remove the unzipped fastq file to save space
    rm ${prefix}.unzipped.fastq

    maln_file=\$(ls ${prefix}.${ref_prefix}.maln.*)
    iteration=\$(echo "\${maln_file}" | grep -oE '[0-9]+\$')

    ## convert format
    ma -M \${maln_file} \\
        -f 41 \\
        > \${maln_file}.41 \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MIA: \$(echo \$(mia 2>&1) | grep -o 'V [0-9.]*')
    END_VERSIONS
    """
}