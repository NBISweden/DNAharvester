process SPLIT_FASTQ {
    tag "$meta.id"
    label 'process_split_fastq'

    conda "conda-forge::gawk=5.3.0"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/gawk:5.3.0--18a2e1510199122d' :
        'community.wave.seqera.io/library/gawk:5.3.0--180f75ae8b0ce739' }"

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
    def split_readlen = params.split_readlen ?: 70

    if (meta.single_end) {
        """
        echo -n | gzip > ${prefix}_short.fastq.gz
        echo -n | gzip > ${prefix}_long.fastq.gz

        zcat ${reads} | awk -v minlen=${split_readlen} '
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
    } else {
        """
        echo -n | gzip > ${prefix}_R1_short.fastq.gz
        echo -n | gzip > ${prefix}_R1_long.fastq.gz
        echo -n | gzip > ${prefix}_R2_short.fastq.gz
        echo -n | gzip > ${prefix}_R2_long.fastq.gz

        ### Process Read 1
        zcat ${reads[0]} | awk -v minlen=${split_readlen} '
        {
            if(NR%4==1) header=\$0;
            else if(NR%4==2) seq=\$0;
            else if(NR%4==3) plus=\$0;
            else if(NR%4==0) {
                qual=\$0;
                if(length(seq) < minlen)
                    print header "\\n" seq "\\n" plus "\\n" qual | "gzip > ${prefix}_R1_short.fastq.gz";
                else
                    print header "\\n" seq "\\n" plus "\\n" qual | "gzip > ${prefix}_R1_long.fastq.gz";
            }
        }'

        ### Process Read 2
        zcat ${reads[1]} | awk -v minlen=${split_readlen} '
        {
            if(NR%4==1) header=\$0;
            else if(NR%4==2) seq=\$0;
            else if(NR%4==3) plus=\$0;
            else if(NR%4==0) {
                qual=\$0;
                if(length(seq) < minlen)
                    print header "\\n" seq "\\n" plus "\\n" qual | "gzip > ${prefix}_R2_short.fastq.gz";
                else
                    print header "\\n" seq "\\n" plus "\\n" qual | "gzip > ${prefix}_R2_long.fastq.gz";
            }
        }'

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            awk: \$(awk --version 2>&1 | head -1 | awk '{print \$1, \$2}')
        END_VERSIONS
        """
    }
}
