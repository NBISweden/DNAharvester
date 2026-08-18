process DEEPVARIANT_CALL {
    tag "$meta.id"
    label 'process_deepvariant'

    container 'google/deepvariant:1.10.0'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta_fasta), path(fasta)
    tuple val(meta_fai), path(fai)

    output:
    tuple val(meta), path("*_raw.vcf.gz")     , emit: raw_vcf
    tuple val(meta), path("*_raw.vcf.gz.tbi") , optional: true, emit: raw_vcf_tbi
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def model_type = task.ext.model_type ?: "${params.deepvariant_model_type}"

    """
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=${model_type} \\
        --ref=${fasta} \\
        --reads=${bam} \\
        --output_vcf=${prefix}_raw.vcf.gz \\
        --num_shards=${task.cpus} \\
        --disable_small_model=true \\
        --vcf_stats_report=true \\
        --logging_dir=logs \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: 1.10.0
    END_VERSIONS
    """
}
