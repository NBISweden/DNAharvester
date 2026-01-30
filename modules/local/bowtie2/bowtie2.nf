process BOWTIE2 {
    tag "$meta.id"
    label 'process_bowtie2'

    conda "bioconda::bowtie2=2.5.4, bioconda::samtools=1.20"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bowtie2_samtools:4f429f1a3a2a870a' :
        'community.wave.seqera.io/library/bowtie2_samtools:ea581343bbe90572' }"

    input:
    tuple val(meta), path(reads)
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
    def bowtie2_params =
        meta.sample_type == 'ancient' ? (params.bowtie2_ancient_params ?: '') :
        meta.sample_type == 'modern' ? (params.bowtie2_modern_params ?: '') : ''
    // Parse read group: @RG\tID:xxx\tSM:xxx\tPL:xxx\tLB:xxx -> --rg-id 'xxx' --rg 'SM:xxx' --rg 'PL:xxx' --rg 'LB:xxx'
    def rg_args = ''
    def parts = meta.read_group.toString().replaceAll('@RG\\\\t', '').tokenize('\\t')
    parts.each { part ->
        if (part.startsWith('ID:')) {
            rg_args = "--rg-id '${part.minus('ID:')}'"
        } else {
            rg_args = "${rg_args} --rg '${part}'"
        }
    }

    if (meta.single_end) {
        // Single-end mapping
        """
        INDEX=`find -L ./ -maxdepth 2 -name "${reference}.1.bt2" | sed 's/\\.1.bt2\$//'`

        bowtie2 \\
            ${args} \\
            ${rg_args} \\
            -p ${task.cpus} \\
            -x \${INDEX} \\
            -U ${reads} | \\
        samtools sort -@ ${task.cpus} -o ${prefix}.${ref_prefix}.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bowtie2: \$(bowtie2 --version 2>&1 | grep -i version | head -n 1 | sed 's/.*version //')
            samtools: \$(samtools --version 2>&1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        // Paired-end mapping
        """
        INDEX=`find -L ./ -maxdepth 2 -name "${reference}.1.bt2" | sed 's/\\.1.bt2\$//'`

        bowtie2 \\
            ${args} \\
            ${rg_args} \\
            -p ${task.cpus} \\
            -x \${INDEX} \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} | \\
        samtools sort -@ ${task.cpus} -o ${prefix}.${ref_prefix}.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bowtie2: \$(bowtie2 --version 2>&1 | grep -i version | head -n 1 | sed 's/.*version //')
            samtools: \$(samtools --version 2>&1 | head -n 1 | sed 's/samtools //')
        END_VERSIONS
        """
    }
}
