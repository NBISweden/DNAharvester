process SAMTOOLS_CHECK_SQ {
    tag "$meta.id"
    label 'process_samplesheet_check'

    conda "bioconda::htslib=1.21 bioconda::samtools=1.21"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9edc2564215d5cd137a8b25ca8a311600987186d406b092022444adf3c4447f7/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.21--6cb89bfd40cbaabf' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path(bam, includeInputs: true)    , emit: bam
    tuple val(meta), path("${meta.id}.sq.tsv")         , emit: sq
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Every reference sequence (name and length) must be present in the BAM @SQ header.
    // Extra @SQ entries (e.g. decoy sequences from competitive mapping) are allowed.
    """
    samtools view -H ${bam} \\
        | awk -F'\\t' '\$1 == "@SQ" { sn = ""; ln = ""; for (i = 2; i <= NF; i++) { if (\$i ~ /^SN:/) sn = substr(\$i, 4); if (\$i ~ /^LN:/) ln = substr(\$i, 4) }; print sn"\\t"ln }' \\
        > ${meta.id}.sq.tsv

    awk -F'\\t' 'NR == FNR { sq[\$1] = \$2; next }
        !(\$1 in sq) { print \$1"\\tnot found in BAM header"; next }
        sq[\$1] != \$2 { print \$1"\\tlength differs (reference: "\$2", BAM: "sq[\$1]")" }' \\
        ${meta.id}.sq.tsv ${fai} > missing_sq.tsv

    if [ -s missing_sq.tsv ]; then
        echo "ERROR: ${bam} was not mapped to the provided reference genome. \$(wc -l < missing_sq.tsv) reference sequence(s) do not match the BAM @SQ header:" >&2
        head -n 20 missing_sq.tsv >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
