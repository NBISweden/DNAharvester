#! /usr/bin/env nextflow

include { FASTQC as FASTQC_RAW         } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_PROCESSED   } from '../../../modules/nf-core/fastqc/main'
include { MAPDAMAGE2                   } from '../../../modules/nf-core/mapdamage2/main'
include { MULTIQC                      } from '../../../modules/nf-core/multiqc/main'

workflow DATA_QC {
    take:
    reference
    raw_reads       // paired-end reads or single-end reads
    processed_reads // merged paired-end reads or trimmed single-end reads
    fastp_json      // read statistic files from FastP
    bam

    main:
    ch_versions                              = Channel.empty()

    FASTQC_RAW ( raw_reads )
    FASTQC_PROCESSED ( processed_reads )
    ch_versions                              = ch_versions.mix(FASTQC_PROCESSED.out.versions)

    MAPDAMAGE2 ( bam, reference )
    ch_versions                              = ch_versions.mix(MAPDAMAGE2.out.versions)

    // Run MultiQC on output
    ch_multiqc_files                         = FASTQC_RAW.out.zip.map{ meta, qcfile -> qcfile }.mix(
                                                FASTQC_PROCESSED.out.zip.map{ meta, qcfile -> qcfile },
                                                fastp_json.map{ meta, fastp_json -> fastp_json },
                                                MAPDAMAGE2.out.pgtoa_freq.map{ meta, pgtoa_freq -> pgtoa_freq },
                                                MAPDAMAGE2.out.pctot_freq.map{ meta, pctot_freq -> pctot_freq },
                                                MAPDAMAGE2.out.lgdistribution.map{ meta, lgdistribution -> lgdistribution },
                                                ).collect()
    ch_multiqc_config                        = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_extra_config                  = params.multiqc_extra_config ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                          = params.multiqc_logo ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(MULTIQC.out.versions)

    emit:
    fastqc_raw_html                          = FASTQC_RAW.out.html                          // channel: [ val(meta), path(html) ]
    fastqc_raw_zip                           = FASTQC_RAW.out.zip                           // channel: [ val(meta), path(zip) ]
    fastqc_processed_html                    = FASTQC_PROCESSED.out.html                    // channel: [ val(meta), path(html) ]
    fastqc_processed_zip                     = FASTQC_PROCESSED.out.zip                     // channel: [ val(meta), path(zip) ]
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
    multiqc_report                           = MULTIQC.out.report.toList()                  // channel: [ val(meta), path(report) ]
    versions                                 = ch_versions                                  // channel: [ versions.yml ]
}