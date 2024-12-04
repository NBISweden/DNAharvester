process AMBER {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::matplotlib=3.8.0 bioconda::pysam=0.21.0"
    container "quay.io/biocontainers/mulled-v2-ecefa487ac2b5705c340ca872c9bb191102f0000:ce2ab207ce57b290cf5fa5816192018e1caa9fb4-0" // mulled image from Biocontainers including the required packages and more

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
