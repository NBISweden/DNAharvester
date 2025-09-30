process ESTIMATE_READ_LEN_CUTOFF {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::kneed=0.8.5"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/kneed:0.8.5--b648f8184207cf07' :
        'community.wave.seqera.io/library/kneed:0.8.5--5e26e0f3d57da086' }"

    input:
    tuple val(meta), path(amber_txt)

    output:
    tuple val(meta), path("*_read_len_cutoff.txt")       , emit: read_len_cutoff
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in {{ name }}/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def curve = task.ext.curve ?: ((params.mapping_tool == 'bwa-aln' || params.mapping_tool == 'bwa-aln-mem') ? 'convex' :
                (params.mapping_tool == 'bowtie2' ? 'concave' : null))
    def direction = task.ext.direction ?: ((params.mapping_tool == 'bwa-aln' || params.mapping_tool == 'bwa-aln-mem') ? 'decreasing' :
                    (params.mapping_tool == 'bowtie2' ? 'increasing' : null))

    """
    ## get the number of reads from filename
    number_reads=\$(cat ${amber_txt} | awk -F': ' '/^sample:/ {print \$2}' | cut -d'_' -f4)

    if [ "\$number_reads" -gt 10000 ]; then
        select_read_len_cutoff_amber.py \\
            ${amber_txt} \\
            ${curve} \\
            ${direction} \\
            > ${prefix}_read_len_cutoff.txt
    else
        echo "Insufficient mapped reads to estimate cutoff (\$number_reads). Using default value of: 30" > ${prefix}_read_len_cutoff.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}