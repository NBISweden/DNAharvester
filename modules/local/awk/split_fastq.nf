process SPLIT_FASTQ {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::gawk=5.3.0"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyp`roject.org/singularity/gawk%3A5.3.1' :
        'docker://stagex/gawk:5.3.0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_short.fastq.gz"), emit: short_reads
    tuple val(meta), path("*_long.fastq.gz") , emit: long_reads
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def minlen = params.split_minlen ?: 35 // Default to 35 if not set in config

    """
    zcat ${reads} | awk -v minlen=${minlen} '
    {
        if(NR%4==1) header=\$0;
        else if(NR%4==2) seq=\$0;
        else if(NR%4==3) plus=\$0;
        else if(NR%4==0) {
            qual=\$0;
            if(length(seq) < minlen)
                print header "\\n" seq "\\n" plus "\\n" qual | "gzip > ${prefix}_short.fastq.gz";
            else
                print header "\\n" seq "\\n" plus "\\n" qual | "gzip > ${prefix}_long.fastq.gz";
        }
    }'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1 | awk '{print \$1, \$2}')
    END_VERSIONS
    """
}
