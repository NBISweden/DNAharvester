#! /usr/bin/env nextflow
process LEEHOM {
    tag "$meta.id"
    label 'process_leehom'

    conda "bioconda::leehom=1.2.15"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/leehom:1.2.15--TODO' :
        'community.wave.seqera.io/library/leehom:1.2.15--TODO' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.se.leehom.fastq.gz')          , optional:true, emit: reads_se
    tuple val(meta), path('*.merged.leehom.fastq.gz')      , optional:true, emit: reads_merged
    tuple val(meta), path('*.unmerged-R*.leehom.fastq.gz') , optional:true, emit: reads_unmerged
    tuple val(meta), path('*.fail.leehom.fastq.gz')        , optional:true, emit: reads_fail
    tuple val(meta), path('*.json')                        , emit: json
    tuple val(meta), path('*.log')                         , emit: log
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def min_length = task.ext.min_length ?: 0
    // Shared helpers: count reads in a fastq(.gz) file, and drop reads shorter
    // than min_length while keeping R1/R2 pairs in sync (leeHom has no built-in
    // minimum-length filter, unlike fastp's -l flag)
    def helper_functions = '''
    count_reads() {
        [ -f "$1" ] && zcat "$1" | awk 'END{print NR/4}' || echo 0
    }

    filter_se_length() {
        local in=$1 out=$2 minlen=$3
        if [ ! -f "$in" ]; then : | gzip -c > "$out"; return; fi
        if [ "$minlen" -le 0 ]; then cp "$in" "$out"; return; fi
        zcat "$in" | paste - - - - | awk -F'\\t' -v len="$minlen" 'length($2) >= len' | tr '\\t' '\\n' | gzip -c > "$out"
    }

    filter_pe_length() {
        local r1_in=$1 r2_in=$2 r1_out=$3 r2_out=$4 minlen=$5
        if [ ! -f "$r1_in" ] || [ ! -f "$r2_in" ]; then
            : | gzip -c > "$r1_out"
            : | gzip -c > "$r2_out"
            return
        fi
        if [ "$minlen" -le 0 ]; then
            cp "$r1_in" "$r1_out"
            cp "$r2_in" "$r2_out"
            return
        fi
        local r1tmp="${r1_out%.gz}.tmp"
        local r2tmp="${r2_out%.gz}.tmp"
        paste <(zcat "$r1_in" | paste - - - -) <(zcat "$r2_in" | paste - - - -) | \\
            awk -F'\\t' -v len="$minlen" -v r1tmp="$r1tmp" -v r2tmp="$r2tmp" '
                length($2) >= len && length($6) >= len {
                    print $1"\\n"$2"\\n"$3"\\n"$4 > r1tmp
                    print $5"\\n"$6"\\n"$7"\\n"$8 > r2tmp
                }
            '
        [ -f "$r1tmp" ] && gzip -c "$r1tmp" > "$r1_out" || (: | gzip -c > "$r1_out")
        [ -f "$r2tmp" ] && gzip -c "$r2tmp" > "$r2_out" || (: | gzip -c > "$r2_out")
        rm -f "$r1tmp" "$r2tmp"
    }
    '''
    if (meta.single_end) {
    """
    $helper_functions

    [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz
    leeHom \\
        -fq1 ${prefix}.fastq.gz \\
        -fqo ${prefix}.leehom \\
        --ancientdna \\
        -t $task.cpus \\
        $args \\
        2> ${prefix}.leehom.log

    filter_se_length ${prefix}.leehom.fq.gz ${prefix}.se.leehom.fastq.gz $min_length
    [ -f ${prefix}.leehom.fail.fq.gz ] && mv ${prefix}.leehom.fail.fq.gz ${prefix}.fail.leehom.fastq.gz

    reads_before=\$(count_reads ${prefix}.fastq.gz)
    reads_after=\$(count_reads ${prefix}.se.leehom.fastq.gz)
    cat <<-END_JSON > ${prefix}.leehom.json
    {"read1_before_filtering":{"total_reads":\${reads_before}},"summary":{"after_filtering":{"total_reads":\${reads_after}}}}
    END_JSON

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leehom: \$(leeHom --version 2>&1)
    END_VERSIONS
    """
    } else {
    """
    $helper_functions

    [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
    leeHom \\
        -fq1 ${prefix}_1.fastq.gz \\
        -fq2 ${prefix}_2.fastq.gz \\
        -fqo ${prefix}.leehom \\
        --ancientdna \\
        -t $task.cpus \\
        $args \\
        2> ${prefix}.leehom.log

    filter_se_length ${prefix}.leehom.fq.gz ${prefix}.merged.leehom.fastq.gz $min_length
    filter_pe_length ${prefix}.leehom_r1.fq.gz ${prefix}.leehom_r2.fq.gz ${prefix}.unmerged-R1.leehom.fastq.gz ${prefix}.unmerged-R2.leehom.fastq.gz $min_length
    [ -f ${prefix}.leehom.fail.fq.gz ]    && mv ${prefix}.leehom.fail.fq.gz    ${prefix}.fail.leehom.fastq.gz
    [ -f ${prefix}.leehom_r1.fail.fq.gz ] && mv ${prefix}.leehom_r1.fail.fq.gz ${prefix}.unmerged-R1.fail.leehom.fastq.gz
    [ -f ${prefix}.leehom_r2.fail.fq.gz ] && mv ${prefix}.leehom_r2.fail.fq.gz ${prefix}.unmerged-R2.fail.leehom.fastq.gz

    reads_before=\$(count_reads ${prefix}_1.fastq.gz)
    reads_after=\$(( \$(count_reads ${prefix}.merged.leehom.fastq.gz) + \$(count_reads ${prefix}.unmerged-R1.leehom.fastq.gz) + \$(count_reads ${prefix}.unmerged-R2.leehom.fastq.gz) ))
    cat <<-END_JSON > ${prefix}.leehom.json
    {"read1_before_filtering":{"total_reads":\${reads_before}},"summary":{"after_filtering":{"total_reads":\${reads_after}}}}
    END_JSON

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leehom: \$(leeHom --version 2>&1)
    END_VERSIONS
    """
    }
}
