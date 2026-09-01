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
    def readlength = (params.readlength == 'auto') ? 20 : params.readlength
    // Custom adapter sequences (optional) - left empty, AdapterRemoval3 falls back to its own defaults/auto-detection
    def adapter1 = params.adapter1 ? "--adapter1 ${params.adapter1}" : ''
    def adapter2 = params.adapter2 ? "--adapter2 ${params.adapter2}" : ''
    // AdapterRemoval3 only merges overlapping mates when explicitly asked to, same as fastp
    def merge_reads_enabled = params.merge_reads.toBoolean()
    def merge_cmd           = merge_reads_enabled ? '--merge' : ''
    def pe_r1_name          = merge_reads_enabled ? 'unmerged-R1' : 'pe.R1'
    def pe_r2_name          = merge_reads_enabled ? 'unmerged-R2' : 'pe.R2'
    // reads_after (below) must mirror exactly what fastq_processing keeps in processed_reads
    // (see subworkflows/local/fastq_processing/main.nf), so the reported read count isn't
    // inflated by categories the pipeline actually discards, e.g. unmerged/singleton reads
    // when keep_unmerged_reads is false
    def keep_unmerged_reads_enabled = params.keep_unmerged_reads.toBoolean()
    def kept_pe_files = []
    if (merge_reads_enabled) {
        kept_pe_files << "${prefix}.merged.adapterremoval.fastq.gz"
        if (keep_unmerged_reads_enabled) {
            kept_pe_files << "${prefix}.${pe_r1_name}.adapterremoval.fastq.gz"
            kept_pe_files << "${prefix}.${pe_r2_name}.adapterremoval.fastq.gz"
            kept_pe_files << "${prefix}.singleton.adapterremoval.fastq.gz"
        }
    } else {
        kept_pe_files << "${prefix}.${pe_r1_name}.adapterremoval.fastq.gz"
        kept_pe_files << "${prefix}.${pe_r2_name}.adapterremoval.fastq.gz"
        if (keep_unmerged_reads_enabled) {
            kept_pe_files << "${prefix}.singleton.adapterremoval.fastq.gz"
        }
    }
    def reads_after_sum = kept_pe_files.collect { "\$(count_reads ${it})" }.join(' + ')
    if (meta.single_end) {
    """
    [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz
    adapterremoval3 \\
        --in-file1 ${prefix}.fastq.gz \\
        --out-prefix ${prefix}.adapterremoval \\
        --threads $task.cpus \\
        --min-length ${readlength} \\
        ${adapter1} \\
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
        --min-length ${readlength} \\
        ${adapter1} \\
        ${adapter2} \\
        $args \\
        2> ${prefix}.adapterremoval.log

    mv ${prefix}.adapterremoval.r1.fastq.gz ${prefix}.${pe_r1_name}.adapterremoval.fastq.gz
    mv ${prefix}.adapterremoval.r2.fastq.gz ${prefix}.${pe_r2_name}.adapterremoval.fastq.gz
    [ -f ${prefix}.adapterremoval.merged.fastq.gz ]    && mv ${prefix}.adapterremoval.merged.fastq.gz    ${prefix}.merged.adapterremoval.fastq.gz
    [ -f ${prefix}.adapterremoval.singleton.fastq.gz ] && mv ${prefix}.adapterremoval.singleton.fastq.gz ${prefix}.singleton.adapterremoval.fastq.gz

    count_reads() { [ -f "\$1" ] && zcat "\$1" | awk 'END{print NR/4}' || echo 0; }
    reads_before=\$(count_reads ${prefix}_1.fastq.gz)
    reads_after=\$(( ${reads_after_sum} ))
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
