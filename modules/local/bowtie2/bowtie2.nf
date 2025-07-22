process BOWTIE2 {
    tag "$meta.id"
    label 'process_bowtie2'

    conda "bioconda::bowtie2=2.5.4, bioconda::samtools=1.20"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bowtie2_samtools:4f429f1a3a2a870a' :
        'quay.io/viascientific/bowtie2-samtools' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '' // ancient DNA parameters are added via args in the configs/modules.config file
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = task.ext.reference ?: "${meta2.id}"
    def ref_prefix = reference.replaceAll(/\.(fasta|fna|fa)$/, '')

    """
    INDEX=`find -L ./ -maxdepth 2 -name "${reference}.1.bt2" | sed 's/\\.1.bt2\$//'`

    bowtie2 --sensitive -p ${task.cpus} -x \${INDEX} -U ${reads} | \
    samtools sort -@ ${task.cpus} -o ${prefix}.${ref_prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | grep -i version | head -n 1 | sed 's/.*version //')
    END_VERSIONS
    """
}