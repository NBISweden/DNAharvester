process BWA_ALN_MEM {
    tag "$meta.id"
    label 'process_bwa_aln_mem'

    conda "bioconda::bwa=0.7.18 bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bwa_samtools:813d9b5fea3890ec' :
        'community.wave.seqera.io/library/bwa_samtools:3938c84206f62975' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("${meta.id}*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '' // ancient DNA parameters are added via args in the configs/modules.config file
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = task.ext.reference ?: "${meta2.id}"
    def ref_prefix = reference.replaceAll(/\.(fasta|fna|fa)$/, '')
    def read_group = task.ext.read_group ?: "${meta.read_group}"

    """
    INDEX=`find -L ./ -maxdepth 2 -name "${reference}.amb" | sed 's/\\.amb\$//'`

    zcat ${reads} | awk -v minlen=90 '
    {
        if(NR%4==1) header=\$0;
        else if(NR%4==2) seq=\$0;
        else if(NR%4==3) plus=\$0;
        else if(NR%4==0) {
            qual=\$0;
            if(length(seq) < minlen)
                print header "\\n" seq "\\n" plus "\\n" qual >> "short_reads.fastq";
            else
                print header "\\n" seq "\\n" plus "\\n" qual >> "long_reads.fastq";
        }
    }'

    bwa aln \\
        -l 16500 -n 0.01 -o 2 \\
        -t ${task.cpus} \\
        -f short_reads.sai \\
        \${INDEX} \\
        short_reads.fastq

    bwa samse \\
        -r "$read_group" \\
        \${INDEX} \\
        short_reads.sai \\
        short_reads.fastq | \\
        samtools sort -@ ${task.cpus} -O bam - > short_reads.bam

    bwa mem \\
        -k 19 -r 2.5 -L 15 \\
        -t ${task.cpus} \\
        -R "$read_group" \\
        \${INDEX} \\
        long_reads.fastq | \\
        samtools sort -@ ${task.cpus} -O bam - > long_reads.bam

    samtools merge -@ ${task.cpus} -o ${prefix}.${ref_prefix}.bam \\
        short_reads.bam long_reads.bam

    ## remove intermediate files
    rm -f short_reads.fastq short_reads.sai short_reads.bam long_reads.fastq long_reads.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}