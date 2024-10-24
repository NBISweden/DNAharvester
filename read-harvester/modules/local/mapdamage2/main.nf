process MAPDAMAGE2 {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::mapdamage2=2.2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mapdamage2%3A2.2.2--pyr43hdfd78af_0' :
        'quay.io/biocontainers/mapdamage2:2.2.2--pyr43hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path(fasta)

    output:
    tuple val(meta), path("${prefix}/Runtime_log.txt")                                    ,emit: runtime_log
    tuple val(meta), path("${prefix}/Fragmisincorporation_plot.pdf"), optional: true      ,emit: fragmisincorporation_plot
    tuple val(meta), path("${prefix}/Length_plot.pdf"), optional: true                    ,emit: length_plot
    tuple val(meta), path("${prefix}/misincorporation.txt"), optional: true               ,emit: misincorporation
    tuple val(meta), path("${prefix}/lgdistribution.txt"), optional: true                 ,emit: lgdistribution
    tuple val(meta), path("${prefix}/dnacomp.txt"), optional: true                        ,emit: dnacomp
    tuple val(meta), path("${prefix}/Stats_out_MCMC_hist.pdf"), optional: true            ,emit: stats_out_mcmc_hist
    tuple val(meta), path("${prefix}/Stats_out_MCMC_iter.csv"), optional: true            ,emit: stats_out_mcmc_iter
    tuple val(meta), path("${prefix}/Stats_out_MCMC_trace.pdf"), optional: true           ,emit: stats_out_mcmc_trace
    tuple val(meta), path("${prefix}/Stats_out_MCMC_iter_summ_stat.csv"), optional: true  ,emit: stats_out_mcmc_iter_summ_stat
    tuple val(meta), path("${prefix}/Stats_out_MCMC_post_pred.pdf"), optional: true       ,emit: stats_out_mcmc_post_pred
    tuple val(meta), path("${prefix}/Stats_out_MCMC_correct_prob.csv"), optional: true    ,emit: stats_out_mcmc_correct_prob
    tuple val(meta), path("${prefix}/dnacomp_genome.csv"), optional: true                 ,emit: dnacomp_genome
    tuple val(meta), path("${prefix}/*rescaled.bam"), optional: true                      ,emit: rescaled
    tuple val(meta), path("${prefix}/5pCtoT_freq.txt"), optional: true                    ,emit: pctot_freq
    tuple val(meta), path("${prefix}/3pGtoA_freq.txt"), optional: true                    ,emit: pgtoa_freq
    tuple val(meta), path("${prefix}/*.fasta"), optional: true                            ,emit: fasta
    tuple val(meta), path("${prefix}/"), optional: true                                   ,emit: folder
    path "versions.yml",emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mapDamage \\
            $args \\
            -d $prefix \\
            -i $bam \\
            -r $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mapdamage2: \$(echo \$(mapDamage --version))
    END_VERSIONS
    """
}
