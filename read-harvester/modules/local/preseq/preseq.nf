process PRESEQ {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::preseq=3.2.0"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/preseq:79160ee386f5eed6' :
        'community.wave.seqera.io/library/preseq:3.2.0--2789d8b704b33613' }"

    input:
    tuple val(meta), path(bam), path(index)

    output:
    tuple val(meta), path("*preseq.txt")    , emit: preseq_txt
    path "versions.yml"                     , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    preseq lc_extrap -B -o ${prefix}_preseq.txt ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preseq: \$(preseq 2>&1 | grep "Version" | sed -e "s/Version: //g")
    END_VERSIONS
    """
}