process GATK_INDEL_REALIGNER {
    tag "$meta.id"
    label 'process_gatk_indel_realigner'

    conda "bioconda::gatk=3.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk%3A3.8--hdfd78af_12' :
        'docker://broadinstitute/gatk3:3.8-1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.realigned.bam")        , emit: bam
    tuple val(meta), path("*.intervals")            , emit: intervals
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? "-Xmx${task.memory.toGiga()}g" : "-Xmx4g"

    """
    echo "Starting GATK Indel Realignment for ${prefix}"
    date

    # Create sequence dictionary if not present
    gatk CreateSequenceDictionary -R ${reference}

    # Step 1: RealignerTargetCreator - identify regions for realignment
    gatk ${avail_mem} \\
        -T RealignerTargetCreator \\
        -R ${reference} \\
        -I ${bam} \\
        -o ${prefix}.intervals \\
        -nt ${task.cpus} \\
        ${args}

    # Step 2: IndelRealigner - perform realignment
    gatk ${avail_mem} \\
        -T IndelRealigner \\
        -R ${reference} \\
        -I ${bam} \\
        -targetIntervals ${prefix}.intervals \\
        -o ${prefix}.realigned.bam

    echo "Finished GATK Indel Realignment for ${prefix}"
    date

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*GATK v//; s/ .*\$//')
    END_VERSIONS
    """
}


