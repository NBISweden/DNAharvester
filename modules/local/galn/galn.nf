// G-ALN: GPU implementation of the BWA-aln algorithm (byte-identical to bwa aln 0.7.19)
// Uses an existing BWA index. Requires an NVIDIA GPU.
process GALN {
    tag "$meta.id"
    label 'process_galn'

    // linux/arm64 image (galn 0.2.7 + samtools 1.25) for the Arrhenius GH200 nodes.
    // Apptainer refuses to pull an arm64 image on an x86 host, so pre-place it in
    // NXF_APPTAINER_CACHEDIR (Nextflow then skips the pull):
    //   bin/fetch_oras_image.sh oras://community.wave.seqera.io/library/galn_samtools:8b100f6efa9d04d1 $NXF_APPTAINER_CACHEDIR
    container 'oras://community.wave.seqera.io/library/galn_samtools:8b100f6efa9d04d1'

    // Expose the GPU inside the container
    containerOptions { workflow.containerEngine in ['apptainer', 'singularity'] ? '--nv' : (workflow.containerEngine == 'docker' ? '--gpus all' : '') }

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
    // galn aln takes the same options as bwa aln
    def aln_params =
        meta.sample_type == 'ancient' ? (params.bwa_aln_ancient_params ?: '') :
        meta.sample_type == 'modern' ? (params.bwa_aln_modern_params ?: '') : ''

    if (meta.single_end) {
        """
        # Locate the BWA index prefix (the container has no `find`)
        for f in ./${reference}.amb ./*/${reference}.amb; do
            [ -e "\$f" ] && INDEX=\${f%.amb} && break
        done
        # galn aln uses all GPUs and CPU cores of the job; samse/sampe follow OMP_NUM_THREADS.
        # Cores come from the scheduler (on Arrhenius: 72 per GPU), so use nproc, not task.cpus.
        export OMP_NUM_THREADS=\$(nproc)

        galn aln \\
            $args \\
            $aln_params \\
            -f ${prefix}.sai \\
            \${INDEX} \\
            ${reads}

        galn samse \\
            -l 0 \\
            $read_group \\
            \${INDEX} \\
            ${prefix}.sai \\
            ${reads} | samtools sort -@ \$(nproc) - > ${prefix}.${ref_prefix}.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            galn: \$(galn 2>&1 | grep -m1 'version' | sed 's/^.*version //')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        # Locate the BWA index prefix (the container has no `find`)
        for f in ./${reference}.amb ./*/${reference}.amb; do
            [ -e "\$f" ] && INDEX=\${f%.amb} && break
        done
        # galn aln uses all GPUs and CPU cores of the job; samse/sampe follow OMP_NUM_THREADS.
        # Cores come from the scheduler (on Arrhenius: 72 per GPU), so use nproc, not task.cpus.
        export OMP_NUM_THREADS=\$(nproc)

        galn aln \\
            $args \\
            $aln_params \\
            -f ${prefix}_1.sai \\
            \${INDEX} \\
            ${reads[0]}

        galn aln \\
            $args \\
            $aln_params \\
            -f ${prefix}_2.sai \\
            \${INDEX} \\
            ${reads[1]}

        galn sampe \\
            -l 0 \\
            $read_group \\
            \${INDEX} \\
            ${prefix}_1.sai ${prefix}_2.sai \\
            ${reads[0]} ${reads[1]} | samtools sort -@ \$(nproc) - > ${prefix}.${ref_prefix}.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            galn: \$(galn 2>&1 | grep -m1 'version' | sed 's/^.*version //')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}
