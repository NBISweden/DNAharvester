process BWA_ALN {
    tag "$meta.id"
    label 'process_bwa_aln'

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
    def read_group = meta.read_group ? "-r '${meta.read_group}'" : ""
    def bwa_aln_params =
        meta.sample_type == 'ancient' ? (params.bwa_aln_ancient_params ?: '') :
        meta.sample_type == 'modern' ? (params.bwa_aln_modern_params ?: '') : ''

    if (meta.single_end) {
        """
        INDEX=`find -L ./ -maxdepth 2 -name "${reference}.amb" | sed 's/\\.amb\$//'`

        bwa aln \\
            $args \\
            $bwa_aln_params \\
            -t $task.cpus \\
            \${INDEX} \\
            ${reads} > ${prefix}.sai

        bwa samse \\
            $read_group \\
            \${INDEX} \\
            ${prefix}.sai \\
            ${reads} | samtools sort -@ ${task.cpus} - > ${prefix}.${ref_prefix}.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        INDEX=`find -L ./ -maxdepth 2 -name "${reference}.amb" | sed 's/\\.amb\$//'`

        bwa aln \\
            $args \\
            $bwa_aln_params \\
            -t $task.cpus \\
            \${INDEX} \\
            ${reads[0]} > ${prefix}_1.sai

        bwa aln \\
            $args \\
            $bwa_aln_params \\
            -t $task.cpus \\
            \${INDEX} \\
            ${reads[1]} > ${prefix}_2.sai

        bwa sampe \\
            $read_group \\
            \${INDEX} \\
            ${prefix}_1.sai ${prefix}_2.sai \\
            ${reads[0]} ${reads[1]} | samtools sort -@ ${task.cpus} - > ${prefix}.${ref_prefix}.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}