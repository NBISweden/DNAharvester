process ANGSD_HAPLOTOPLINK {
    tag "$meta.id"
    label 'process_angsd_haplotoplink'

    conda "bioconda::angsd=0.939"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/angsd:0.939--h468462d_0':
        'quay.io/biocontainers/angsd:0.939--h468462d_0' }"

    input:
    tuple val(meta), path(haplo)

    output:
    tuple val(meta), path("*.tfam")                    , emit: tfam
    tuple val(meta), path("*.tped")                    , emit: tped
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    haploToPlink \\
        ${prefix}.haplo.gz \\
        ${prefix} &&

    # Modify sample name in *.tfam output
    sed -i 's/ind0/${prefix}/g' ${prefix}.tfam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        angsd: \$(echo \$(angsd 2>&1) | grep version | head -n 1 | sed 's/.*version: //g;s/ .*//g')
    END_VERSIONS
    """
}