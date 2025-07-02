process RM_TRANSITIONS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::pysam=0.23.3 conda-forge::tqdm=4.67.1"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_tqdm:64064577e6a6659c' :
        'oras://community.wave.seqera.io/library/pysam_tqdm:64064577e6a6659c' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)

    output:
    tuple val(meta), path("*rm_trans.bam")   , emit: rm_trans_bam
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def library_type = meta.library_type


    adna_sslib_damage_removal = "https://raw.githubusercontent.com/bilalbioinfo/aDNA_Damage/refs/heads/main/adna_sslib_damage_removal.py"
    """
    if [ ! -f ${projectDir}/bin/adna_sslib_damage_removal.py ]; then
        wget $adna_sslib_damage_removal &&
        chmod +x adna_sslib_damage_removal.py &&
        mv adna_sslib_damage_removal.py ${projectDir}/bin/
    fi

    adna_sslib_damage_removal.py \\
        --bamfile $bam \\
        --output ${prefix}.rm_trans.bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}