process SAMTOOLS_VIEW_SUBSAMPLE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.19.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)
    path(fasta)

    output:
    tuple val(meta), path("${prefix}.${subsample}_reads.bam") , optional:true, emit: bam
    tuple val(meta), path("${prefix}.${subsample}_reads.cram"), optional:true, emit: cram
    path "versions.yml",                                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def subsample = task.ext.subsample ?: '1'
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    def file_type = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt bam") ? "bam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    input.getExtension()
    """
    total=\$(samtools view -c $bam) # total number of reads, not pairs (may include unmapped and duplicated multi-aligned reads)
    frac=\$(awk -v s=$subsample -v t=\$total "BEGIN {print s/t}") # fraction of reads to keep

    samtools \\
        view \\
        --threads ${task.cpus-1} \\
        ${reference} \\
        -s \$frac \\
        $args \\
        -o ${prefix}.${subsample}_reads.${file_type} \\
        $bam \\
        $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt bam") ? "bam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    input.getExtension()
    if ("$input" == "${prefix}.${file_type}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    def index = args.contains("--write-index") ? "touch ${prefix}.csi" : ""

    """
    touch ${prefix}.${file_type}
    ${index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}