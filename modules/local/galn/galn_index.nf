// G-ALN index: builds the BWA index on the GPU (byte-identical to `bwa index`, same files).
// Drop-in replacement for BWA_INDEX: same inputs/outputs, and storeDir skips it if the index exists.
// The builder can be chosen with ext.args, e.g. '-a gpu-lowmem' or '-a cpu' (default: auto).
// Uses all GPUs and CPU cores of the job.
process GALN_INDEX {
    tag "$fasta"
    label 'process_galn_index'

    // linux/arm64 image (galn 0.2.7 + samtools 1.25), see modules/local/galn/galn.nf
    container 'oras://community.wave.seqera.io/library/galn_samtools:8b100f6efa9d04d1'

    // Expose the GPU inside the container
    containerOptions { workflow.containerEngine in ['apptainer', 'singularity'] ? '--nv' : (workflow.containerEngine == 'docker' ? '--gpus all' : '') }

    storeDir "${output_dir}"

    input:
    tuple val(meta2), path(fasta)
    val (output_dir)

    output:
    tuple val(meta2), path("${meta2.id}.{amb,ann,bwt,pac,sa}", arity: '5')  , emit: index
    tuple val(meta2), val(output_dir)                                       , emit: index_dir
    path "versions.yml"                                                     , emit: versions, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    galn \\
        index \\
        $args \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        galn: \$(galn 2>&1 | grep -m1 'version' | sed 's/^.*version //')
    END_VERSIONS
    """
}
