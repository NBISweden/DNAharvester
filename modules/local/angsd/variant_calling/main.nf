process ANGSD_VARIANT_CALLING {
    tag "$meta.id"
    label 'process_angsd_variant_calling'

    conda "bioconda::angsd=0.939"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/angsd:0.939--h468462d_0':
        'quay.io/biocontainers/angsd:0.939--h468462d_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(fai)

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


    """
    ls -1 *.bam > ${prefix}.bamlist.txt

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