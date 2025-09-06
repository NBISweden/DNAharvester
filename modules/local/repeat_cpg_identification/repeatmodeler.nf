process REPEATMODELER {
    tag "$fasta"
    label 'process_repeatmodeler'

    conda "bioconda::repeatmodeler=2.0.6"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/repeatmodeler:2.0.6--33ecba32dc152693' :
        'community.wave.seqera.io/library/repeatmodeler:2.0.6--64a830a44f180fb9' }"

    input:
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta2), path ("*.upper.fasta")                    , emit: upper_fasta
    tuple val(meta2), path ("RM_*.*/consensi.fa.classified")    , emit: consensi
    tuple val(meta2), path ("RM_*.*/families-classified.stk")   , emit: families
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta2.id}"
    def args = task.ext.args ?: ''

    """
    ## changing the reference fasta to upper case
    awk '{{ if (\$0 !~ />/) {{print toupper(\$0)}} else {{print \$0}} }}' ${fasta} > ${prefix}.upper.fasta

    ## build repeat database
    BuildDatabase -name ${prefix}_database ${prefix}.upper.fasta

    ## run RepeatModeler
    RepeatModeler -engine ncbi -threads ${task.cpus} -database ${prefix}_database -quick

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RepeatModeler: \$(echo \$(RepeatModeler --version 2>&1) | awk '{print \$3}')
    END_VERSIONS
    """
}


