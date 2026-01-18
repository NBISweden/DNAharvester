process GATK_INDEL_REALIGNER {
    tag "$meta.id"
    label 'process_gatk_indel_realigner'

    conda "bioconda::gatk=3.8"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'docker://broadinstitute/gatk3:3.8-1' :
        'docker://broadinstitute/gatk3:3.8-1' }" // same containers for both all Engines since jar is located at difference places

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(reference)
    tuple val(meta2), path(fai)
    tuple val(meta3), path(dict)

    output:
    tuple val(meta), path("*.realigned.bam")        , emit: realigned_bam
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? (task.memory.toGiga()).toInteger() : 4
    def ref_prefix = task.ext.ref_prefix ?: "${meta2.id}".replaceAll(/\.(fasta|fna|fa)$/, '')

    """
    ### Step 1: RealignerTargetCreator - identify regions for realignment
    java -Xmx${avail_mem}g -jar /usr/GenomeAnalysisTK.jar \\
        -T RealignerTargetCreator \\
        -R ${reference} \\
        -I ${bam} \\
        -o ${prefix}.intervals \\
        ${args}

    ### Step 2: IndelRealigner
    java -Xmx${avail_mem}g -jar /usr/GenomeAnalysisTK.jar \\
        -T IndelRealigner \\
        -R ${reference} \\
        -I ${bam} \\
        -targetIntervals ${prefix}.intervals \\
        -o ${prefix}.${ref_prefix}.realigned.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: 3.8
    END_VERSIONS
    """
}


