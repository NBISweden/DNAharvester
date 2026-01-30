process BWA_MEM {
    tag "$meta.id"
    label 'process_bwa_mem'

    conda "bioconda::bwa=0.7.18 bioconda::samtools=1.20"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bwa_samtools:813d9b5fea3890ec' :
        'community.wave.seqera.io/library/bwa_samtools:3938c84206f62975' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = task.ext.reference ?: "${meta2.id}"
    def ref_prefix = reference.replaceAll(/\.(fasta|fna|fa)$/, '')
    def read_group = meta.read_group ? "${meta.read_group}" : ""
    def bwa_mem_params =
        meta.sample_type == 'ancient' ? (params.bwa_mem_ancient_params ?: '') :
        meta.sample_type == 'modern' ? (params.bwa_mem_modern_params ?: '') : ''

    if (meta.single_end) {
        """
        INDEX=`find -L ./ -maxdepth 2 -name "${reference}.amb" | sed 's/\\.amb\$//'`

        bwa mem \\
            $args \\
            $bwa_mem_params \\
            -t ${task.cpus} \\
            -R "${read_group}" \\
            \${INDEX} \\
            ${reads} | \\
        samtools sort -@ ${task.cpus} -o ${prefix}.${ref_prefix}.bam -

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        INDEX=`find -L ./ -maxdepth 2 -name "${reference}.amb" | sed 's/\\.amb\$//'`

        bwa mem \\
            $args \\
            $bwa_mem_params \\
            -t ${task.cpus} \\
            -R "${read_group}" \\
            \${INDEX} \\
            ${reads[0]} ${reads[1]} | \\
        samtools sort -@ ${task.cpus} -o ${prefix}.${ref_prefix}.bam -

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}
