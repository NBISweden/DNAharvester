process BWA_INDEX {
    tag "$fasta"
    label 'process_medium'

    conda "bioconda::bwa=0.7.18"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bwa:0.7.18--4543b4091f454101' :
        'community.wave.seqera.io/library/bwa:0.7.18--324359fbc6e00dba' }"

    storeDir "${params.reference.substring(0, params.reference.lastIndexOf('/'))}"

    input:
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta2), path("*.amb"), path("*.ann"), path("*.bwt"), path("*.pac"), path("*.sa")  , emit: index
    path "versions.yml"                                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bwa \\
        index \\
        $args \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}