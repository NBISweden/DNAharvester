process AMBER {
    tag "$meta.id"
    label 'process_single'

    container "ghcr.io/NBISweden/LTS-L_Dalen_2302/AMBER"

    input:
    tuple val(meta), path(bamfiles)


    output:
    tuple val(meta), path("amber_plot.pdf")                                    , emit: plot

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python AMBER \\
            $args \\
            --bamfiles $bamfiles
    """
}
