#! /usr/bin/env nextflow
process ADAPTERREMOVAL {
    tag "$meta.id"
    label 'process_adapterremoval'

    conda "bioconda::adapterremoval=3.0.1"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/adapterremoval%3A3.0.1--pl5321h0f5e619_0' :
        'oras://community.wave.seqera.io/library/adapterremoval:3.0.1--25445162faa2b6b8' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.se.adapterremoval.fastq.gz')          , optional:true, emit: reads_se
    tuple val(meta), path('*.merged.adapterremoval.fastq.gz')      , optional:true, emit: reads_merged
    tuple val(meta), path('*.pe.R*.adapterremoval.fastq.gz')       , optional:true, emit: reads_pe
    tuple val(meta), path('*.unmerged-R*.adapterremoval.fastq.gz') , optional:true, emit: reads_unmerged
    tuple val(meta), path('*.singleton.adapterremoval.fastq.gz')   , optional:true, emit: reads_singleton
    tuple val(meta), path('*.json')                                , emit: json
    tuple val(meta), path('*.html')                                , optional:true, emit: html
    tuple val(meta), path('*.log')                                 , emit: log
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def prefix = task.ext.prefix ?: "${meta.id}"
    // AdapterRemoval3 only merges overlapping mates when explicitly asked to, same as fastp
    def merge_reads_enabled = params.merge_reads.toBoolean()
    def merge_cmd           = merge_reads_enabled ? '--merge' : ''
    def pe_r1_name          = merge_reads_enabled ? 'unmerged-R1' : 'pe.R1'
    def pe_r2_name          = merge_reads_enabled ? 'unmerged-R2' : 'pe.R2'
    if (meta.single_end) {
    """
    [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz
    adapterremoval3 \\
        --in-file1 ${prefix}.fastq.gz \\
        --out-prefix ${prefix}.adapterremoval \\
        --threads $task.cpus \\
        $args \\
        2> ${prefix}.adapterremoval.log

    mv ${prefix}.adapterremoval.r1.fastq.gz ${prefix}.se.adapterremoval.fastq.gz

    count_reads() { [ -f "\$1" ] && zcat "\$1" | awk 'END{print NR/4}' || echo 0; }
    reads_before=\$(count_reads ${prefix}.fastq.gz)
    reads_after=\$(count_reads ${prefix}.se.adapterremoval.fastq.gz)
    cat <<-END_JSON > ${prefix}.adapterremoval.json
    {"read1_before_filtering":{"total_reads":\${reads_before}},"summary":{"after_filtering":{"total_reads":\${reads_after}}}}
    END_JSON

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        adapterremoval: \$(adapterremoval3 --version 2>&1)
    END_VERSIONS
    """
    } else {
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
    adapterremoval3 \\
        --in-file1 ${prefix}_1.fastq.gz \\
        --in-file2 ${prefix}_2.fastq.gz \\
        --out-prefix ${prefix}.adapterremoval \\
        --threads $task.cpus \\
        $merge_cmd \\
        $args \\
        2> ${prefix}.adapterremoval.log

    mv ${prefix}.adapterremoval.r1.fastq.gz ${prefix}.${pe_r1_name}.adapterremoval.fastq.gz
    mv ${prefix}.adapterremoval.r2.fastq.gz ${prefix}.${pe_r2_name}.adapterremoval.fastq.gz
    [ -f ${prefix}.adapterremoval.merged.fastq.gz ]    && mv ${prefix}.adapterremoval.merged.fastq.gz    ${prefix}.merged.adapterremoval.fastq.gz
    [ -f ${prefix}.adapterremoval.singleton.fastq.gz ] && mv ${prefix}.adapterremoval.singleton.fastq.gz ${prefix}.singleton.adapterremoval.fastq.gz

    count_reads() { [ -f "\$1" ] && zcat "\$1" | awk 'END{print NR/4}' || echo 0; }
    reads_before=\$(count_reads ${prefix}_1.fastq.gz)
    reads_after=\$(( \$(count_reads ${prefix}.${pe_r1_name}.adapterremoval.fastq.gz) + \$(count_reads ${prefix}.${pe_r2_name}.adapterremoval.fastq.gz) + \$(count_reads ${prefix}.merged.adapterremoval.fastq.gz) + \$(count_reads ${prefix}.singleton.adapterremoval.fastq.gz) ))
    cat <<-END_JSON > ${prefix}.adapterremoval.json
    {"read1_before_filtering":{"total_reads":\${reads_before}},"summary":{"after_filtering":{"total_reads":\${reads_after}}}}
    END_JSON

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        adapterremoval: \$(adapterremoval3 --version 2>&1)
    END_VERSIONS
    """
    }
}
