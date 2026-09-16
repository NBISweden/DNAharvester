process CREATE_SEQUENCE_DICTIONARY {
    tag "$meta2.id"
    label 'process_create_sequence_dictionary'

    conda "bioconda::picard"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'docker://broadinstitute/picard:3.4.0' :
        'broadinstitute/picard:3.4.0' }" // same containers for both all Engines since jar is located at difference places

    input:
    tuple val(meta2), path(reference)

    output:
    tuple val(meta2), path("*.dict")                , emit: dict
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Reserve 1GB of headroom below task.memory for the JVM's own overhead (metaspace, thread
    // stacks, GC, JIT code cache, mmap'd CDS archive) - without it, -Xmx == the container's
    // cgroup memory limit and the JVM's real footprint breaches that limit almost immediately
    def avail_mem = task.memory ? Math.max((task.memory.toGiga() - 1).toInteger(), 1) : 4

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

