process RM_TRANSITIONS {
    tag "$meta.id"
    label 'process_rm_transitions'

    conda "bioconda::pysam=0.23.3 conda-forge::tqdm=4.67.1 wget=1.21.4"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_tqdm_wget:c28479ff2635fb54' :
        'community.wave.seqera.io/library/pysam_tqdm_wget:854394b0b7ab4567' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*rm_trans.bam")   , emit: rm_trans_bam
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def library_type = meta.library_type


    adna_sslib_damage_removal = "https://raw.githubusercontent.com/bilalbioinfo/aDNA_Damage/refs/heads/main/deamstrip.py"
    """
    if [ ! -f ${projectDir}/bin/deamstrip.py ]; then
        wget $adna_sslib_damage_removal &&
        chmod +x deamstrip.py &&
        mv deamstrip.py ${projectDir}/bin/
    fi

    deamstrip.py \\
        --bamfile $bam \\
        --library $library_type \\
        --output ${prefix}.rm_trans.bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}