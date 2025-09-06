process BWA_SAMPE {
    tag "$meta.id"
    label 'process_bwa_samse'

    conda "bioconda::bwa=0.7.18 bioconda::samtools=1.20"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bwa_samtools:813d9b5fea3890ec' :
        'community.wave.seqera.io/library/bwa_samtools:3938c84206f62975' }"

    input:
    tuple val(meta), path(read_1), path(read_2), path(sai_1), path(sai_2)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_group = meta.read_group ? "-r '${meta.read_group}'" : ""
    def reference = task.ext.reference ?: "${meta2.id}"
    def ref_prefix = reference.replaceAll(/\.(fasta|fna|fa)$/, '')

    """
    INDEX=`find -L ./ -maxdepth 2 -name "${reference}.amb" | sed 's/\\.amb\$//'`

    bwa sampe \\
        $args \\
        $read_group \\
        \${INDEX} \\
        $sai_1 $sai_2 \\
        $read_1 $read_2 | samtools sort -@ ${task.cpus - 1} -O bam - > ${prefix}.${ref_prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}