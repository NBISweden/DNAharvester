process CREATE_SEQUENCE_DICTIONARY {
    tag "$meta2.id"
    label 'process_create_sequence_dictionary'

    conda "bioconda::picard"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://broadinstitute/picard:3.4.0' :
        'docker://broadinstitute/picard:3.4.0' }" // same container since jar is located at difference places

    input:
    tuple val(meta2), path(reference)

    output:
    tuple val(meta2), path("*.dict")                , emit: dict
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = task.memory ? (task.memory.toGiga()).toInteger() : 4

    """
    # Create sequence dictionary if not present
    java -Xmx${avail_mem}g -jar /usr/picard/picard.jar \\
        CreateSequenceDictionary \\
        REFERENCE=$reference \\
        OUTPUT=${reference.baseName}.dict \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 3.4.0
    END_VERSIONS
    """
}

