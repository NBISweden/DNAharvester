process ANGSD_DOHAPLOCALL {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::angsd=0.939"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/angsd:0.939--h468462d_0':
        'biocontainers/angsd:0.939--h468462d_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.haplo.gz")                , emit: haplo
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ls -1 *.bam > bamlist.txt

    angsd \\
        -bam bamlist.txt \\
        -nThreads ${task.cpus} \\
        -dohaplocall 1 \\
        -doCounts 1 \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        angsd: \$(echo \$(angsd 2>&1) | grep version | head -n 1 | sed 's/.*version: //g;s/ .*//g')
    END_VERSIONS
    """
}