process AMBER {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::matplotlib=3.9.3 bioconda::pysam=0.22.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_matplotlib:9a83292e6b804598' :
        'community.wave.seqera.io/library/pysam_matplotlib:58a92b14d0d8ded9' }"

    input:
    tuple val(meta), path(bam), path(tsv)

    output:
    tuple val(meta), path("*.amber_plot.pdf"), emit: plot
    tuple val(meta), path("*.amber_plot.txt"), emit: txt
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is available at https://github.com/tvandervalk/AMBER/AMBER and has to be placed in bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    AMBER \\
        $args \\
        --bamfiles $tsv \\
        --output ${prefix}.amber_plot

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
