process RM_TRANSITIONS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::matplotlib=3.9.3 bioconda::pysam=0.22.1 conda-forge::wget=1.21.4"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_matplotlib_wget:29119c93f69dc707' :
        'community.wave.seqera.io/library/pysam_matplotlib_wget:a20bf1a7f1b8bebe' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*rm_trans.bam")   , emit: rm_trans_bam
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def library_type = meta.library_type


    adna_sslib_damage_removal = "https://github.com/CpgSthlm/aDNA_Damage/blob/main/adna_sslib_damage_removal.py"
    """
    if [ ! -f ${projectDir}/bin/adna_sslib_damage_removal.py ]; then
        wget $adna_sslib_damage_removal &&
        chmod +x adna_sslib_damage_removal.py &&
        mv adna_sslib_damage_removal.py ${projectDir}/bin/
    fi

    python3
        adna_sslib_damage_removal.py \\
        --bamfiles $bam \\
        --output ${prefix}.rm_trans.bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}