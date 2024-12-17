#! /usr/bin/env nextflow

include { SAMTOOLS_FLAGSTAT          } from '../../../modules/nf-core/samtools/flagstat/main'
include { MAPDAMAGE2                 } from '../../../modules/local/mapdamage2/main'
include { SAMTOOLS_VIEW_SUBSAMPLE    } from '../../../modules/local/samtools/view_subsample/main'
include { CREATE_AMBER_SAMPLESHEET   } from '../../../modules/local/amber/create_amber_samplesheet'
include { AMBER                      } from '../../../modules/local/amber/amber'
include { MULTIQC as MULTIQC_BAM     } from '../../../modules/nf-core/multiqc/main'


workflow RAW_BAM_QC {
    take:
    competitive_reference
    reference
    bam             // bam file from mapping subworkflow
    bai             // bam index file

    main:
    ch_versions                              = Channel.empty()

    ch_samtools_flagstat                     = bam.join(bai)

    SAMTOOLS_FLAGSTAT ( ch_samtools_flagstat )
    ch_versions                              = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    MAPDAMAGE2 ( bam, competitive_reference )
    ch_versions                              = ch_versions.mix(MAPDAMAGE2.out.versions)

    // Run MultiQC on MapDamage and samtools flagstat output
    ch_multiqc_bam_files                     = MAPDAMAGE2.out.folder.map{ meta, folder -> folder }.mix(
                                                SAMTOOLS_FLAGSTAT.out.flagstat.map{ meta, flagstat -> flagstat }).collect()
    ch_multiqc_config                        = params.multiqc_config       ? Channel.fromPath( params.multiqc_config,       checkIfExists: true ) : Channel.empty()
    ch_multiqc_extra_config                  = params.multiqc_extra_config ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                          = params.multiqc_logo         ? Channel.fromPath( params.multiqc_logo,         checkIfExists: true ) : Channel.empty()

    MULTIQC_BAM (
        ch_multiqc_bam_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(MULTIQC_BAM.out.versions)

    SAMTOOLS_VIEW_SUBSAMPLE ( bam, reference )
    ch_versions                              = ch_versions.mix(SAMTOOLS_VIEW_SUBSAMPLE.out.versions)

    CREATE_AMBER_SAMPLESHEET ( SAMTOOLS_VIEW_SUBSAMPLE.out.subsampled_bam )
    ch_versions                              = ch_versions.mix(CREATE_AMBER_SAMPLESHEET.out.versions)

    AMBER (
        SAMTOOLS_VIEW_SUBSAMPLE.out.subsampled_bam.join( CREATE_AMBER_SAMPLESHEET.out.tsv )
    )
    ch_versions                              = ch_versions.mix(AMBER.out.versions)

    emit:
    mapdamage2_fragmisincorporation_plot     = MAPDAMAGE2.out.fragmisincorporation_plot     // channel: [ val(meta), path(fragmisincorporation_plot) ]
    mapdamage2_length_plot                   = MAPDAMAGE2.out.length_plot                   // channel: [ val(meta), path(length_plot) ]
    mapdamage2_misincorporation              = MAPDAMAGE2.out.misincorporation              // channel: [ val(meta), path(misincorporation) ]
    mapdamage2_lgdistribution                = MAPDAMAGE2.out.lgdistribution                // channel: [ val(meta), path(lgdistribution) ]
    mapdamage2_dnacomp                       = MAPDAMAGE2.out.dnacomp                       // channel: [ val(meta), path(dnacomp) ]
    mapdamage2_stats_out_mcmc_hist           = MAPDAMAGE2.out.stats_out_mcmc_hist           // channel: [ val(meta), path(stats_out_mcmc_hist) ]
    mapdamage2_stats_out_mcmc_iter           = MAPDAMAGE2.out.stats_out_mcmc_iter           // channel: [ val(meta), path(stats_out_mcmc_iter) ]
    mapdamage2_stats_out_mcmc_trace          = MAPDAMAGE2.out.stats_out_mcmc_trace          // channel: [ val(meta), path(stats_out_mcmc_trace) ]
    mapdamage2_stats_out_mcmc_iter_summ_stat = MAPDAMAGE2.out.stats_out_mcmc_iter_summ_stat // channel: [ val(meta), path(stats_out_mcmc_iter_summ_stat) ]
    mapdamage2_stats_out_mcmc_post_pred      = MAPDAMAGE2.out.stats_out_mcmc_post_pred      // channel: [ val(meta), path(stats_out_mcmc_post_pred) ]
    mapdamage2_stats_out_mcmc_correct_prob   = MAPDAMAGE2.out.stats_out_mcmc_correct_prob   // channel: [ val(meta), path(stats_out_mcmc_correct_prob) ]
    mapdamage2_dnacomp_genome                = MAPDAMAGE2.out.dnacomp_genome                // channel: [ val(meta), path(dnacomp_genome) ]
    mapdamage2_pctot_freq                    = MAPDAMAGE2.out.pctot_freq                    // channel: [ val(meta), path(pctot_freq) ]
    mapdamage2_pgtoa_freq                    = MAPDAMAGE2.out.pgtoa_freq                    // channel: [ val(meta), path(pgtoa_freq) ]
    mapdamage2_folder                        = MAPDAMAGE2.out.folder                        // channel: [ val(meta), path(folder) ]
    flagstat                                 = SAMTOOLS_FLAGSTAT.out.flagstat               // channel: [ val(meta), path(flagstat) ]
    subsampled_bam                           = SAMTOOLS_VIEW_SUBSAMPLE.out.subsampled_bam   // channel: [ val(meta), path(bam) ]
    amber_plot                               = AMBER.out.plot                               // channel: [ val(meta), path(plot) ]
    amber_txt                                = AMBER.out.txt                                // channel: [ val(meta), path(txt) ]
    multiqc_bam_report                       = MULTIQC_BAM.out.report.toList()              // channel: [ val(meta), path(report) ]
    versions                                 = ch_versions                                  // channel: [ versions.yml ]
}