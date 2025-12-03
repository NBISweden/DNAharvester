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
    tuple val(meta), path("*.sai"), emit: sai
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = task.ext.reference ?: "${meta2.id}"
    def ref_prefix = reference.replaceAll(/\.(fasta|fna|fa)$/, '')
    def bwa_aln_params =
        meta.sample_type == 'ancient' ? (params.bwa_aln_ancient_params ?: '') :
        meta.sample_type == 'modern' ? (params.bwa_aln_modern_params ?: '') : ''

    """
    INDEX=`find -L ./ -maxdepth 2 -name "${reference}.amb" | sed 's/\\.amb\$//'`

    bwa aln \\
        $args \\
        $bwa_aln_params \\
        -t $task.cpus \\
        -f ${prefix}.${ref_prefix}.sai \\
        \${INDEX} \\
        ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}