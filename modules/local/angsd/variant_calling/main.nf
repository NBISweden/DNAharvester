process ANGSD_VARIANT_CALLING {
    tag "$meta.id"
    label 'process_angsd_variant_calling'

    conda "bioconda::angsd=0.939"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/angsd:0.939--h468462d_0':
        'quay.io/biocontainers/angsd:0.939--h468462d_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(fai)
    tuple val(meta2), path(bed_file)

    output:
    tuple val(meta), path("*.arg")                   , emit: angsd_log
    tuple val(meta), path("*.bamlist.txt")           , emit: angsd_bamlist
    tuple val(meta), path("*.geno.gz")               , emit: angsd_geno
    tuple val(meta), path("*.mafs.gz")               , emit: angsd_mafs
    tuple val(meta), path("*.beagle.gz")             , emit: angsd_beagle
    tuple val(meta), path("*.bcf")                   , emit: angsd_bcf
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def filters = task.ext.filters ?: "${params.angsd_filters}"
    def has_regions = !(meta2.id == 'null' && bed_file.name == 'null')
    def regions_file = has_regions ? '-sites "$angsd_sites"' : ""


    """
    ls -1 *.bam > ${prefix}.bamlist.txt

    ${has_regions ? "angsd_sites=\"\$(basename ${bed_file} .bed).angsd\"" : ""}

    # Create regions file for ANGSD if BED file is provided
    # or if repeat-masked BED file is generated via
    # repeat_cpg_identification
    ${has_regions ? "awk '{print \$1\"\\t\"\$2+1\"\\t\"\$3}' ${bed_file} > \"\$angsd_sites\"" : ""}
    ${has_regions ? "angsd sites index \"\$angsd_sites\"" : ""}

    angsd -bam ${prefix}.bamlist.txt \\
        -ref $fasta \\
        -fai $fai \\
        -doMaf 1 \\
        -doMajorMinor 1 \\
        -dogeno 1 \\
        -docounts 1 \\
        -doGlf 2 \\
        -GL 1 \\
        -doPost 1 \\
        -doBcf 1 \\
        ${regions_file} \\
        -nThreads ${task.cpus} \\
        -out ${prefix} \\
        ${filters} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        angsd: \$(echo \$(angsd 2>&1) | grep version | head -n 1 | sed 's/.*version: //g;s/ .*//g')
    END_VERSIONS
    """
}