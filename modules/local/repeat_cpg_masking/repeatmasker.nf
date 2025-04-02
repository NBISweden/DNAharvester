process REPEATMASKER {
    tag "$fasta"
    label 'process_repeatmasker'

    conda "bioconda::repeatmasker=4.1.8"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.8--pl5321hdfd78af_0' :
        'quay.io/biocontainers/repeatmasker:4.1.8--pl5321hdfd78af_0' }"

    input:
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(upper_ref)
    tuple val(meta2), path(consensi)

    output:
    tuple val(meta2), path ("*.out")       , emit: repeatmasker_out
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta2.id}"
    def args = task.ext.args ?: ''
    """
    RepeatMasker -pa ${task.cpus} -xsmall -gccalc \\
        -dir ./ -lib ${consensi} ${upper_ref}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RepeatMasker: \$(echo \$(RepeatMasker -v 2>&1) | awk '{print \$3}')
    END_VERSIONS
    """
}
