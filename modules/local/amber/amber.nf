process AMBER {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::matplotlib=3.9.3 bioconda::pysam=0.22.1 conda-forge::wget=1.21.4"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_matplotlib_wget:29119c93f69dc707' :
        'community.wave.seqera.io/library/pysam_matplotlib_wget:a20bf1a7f1b8bebe' }"

    input:
    tuple val(meta), path(bam), path(tsv)

    output:
    tuple val(meta), path("*.amber.pdf"), optional: true    , emit: plot
    tuple val(meta), path("*.amber.txt")                    , emit: txt
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    amber = "https://raw.githubusercontent.com/tvandervalk/AMBER/refs/heads/main/AMBER"
    """
    if [ ! -f ${projectDir}/bin/AMBER ]; then
        wget $amber &&
        chmod +x AMBER &&
        mv AMBER ${projectDir}/bin/
    fi

    AMBER \\
        $args \\
        --bamfiles $tsv \\
        --output ${prefix}.amber \\
        || echo "No AMBER output was generated, likely due to an insufficient number of mapped reads.\n \\
        Please check the sequencing statistics!." >> ${prefix}.amber_plot.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
