#! /usr/bin/env nextflow

include { SAMTOOLS_FLAGSTAT     as PBQC_FLAGSTAT_MQ_FILTERED_BAM        } from '../../../modules/local/samtools/samtools_flagstat.nf'
include { MULTIQC               as PBQC_MULTIQC_MQ_FILTERED_BAM         } from '../../../modules/nf-core/multiqc/main'
include { SAMTOOLS_FLAGSTAT     as PBQC_FLAGSTAT_RM_SHORT_READS_BAM     } from '../../../modules/local/samtools/samtools_flagstat.nf'
include { MULTIQC               as PBQC_MULTIQC_RM_SHORT_READS_BAM      } from '../../../modules/nf-core/multiqc/main'
include { SAMTOOLS_FLAGSTAT     as PBQC_FLAGSTAT_MERGED_BAM_LIB         } from '../../../modules/local/samtools/samtools_flagstat.nf'
include { MULTIQC               as PBQC_MULTIQC_MERGED_BAM_LIB          } from '../../../modules/nf-core/multiqc/main'
include { PRESEQ                as PBQC_PRESEQ                          } from '../../../modules/local/preseq/preseq'
include { PLOT_PRESEQ           as PBQC_PLOT_PRESEQ                     } from '../../../modules/local/preseq/plot_preseq'
include { SAMTOOLS_FLAGSTAT     as PBQC_FLAGSTAT_DEDUP_LIB              } from '../../../modules/local/samtools/samtools_flagstat.nf'
include { MULTIQC               as PBQC_MULTIQC_DEDUP_LIB               } from '../../../modules/nf-core/multiqc/main'
include { SAMTOOLS_FLAGSTAT     as PBQC_FLAGSTAT_MERGED_BAM_SAMPLE      } from '../../../modules/local/samtools/samtools_flagstat.nf'
include { MULTIQC               as PBQC_MULTIQC_MERGED_BAM_SAMPLE       } from '../../../modules/nf-core/multiqc/main'
include { SAMTOOLS_FLAGSTAT     as PBQC_FLAGSTAT_DEDUP_SAMPLE           } from '../../../modules/local/samtools/samtools_flagstat.nf'
include { MULTIQC               as PBQC_MULTIQC_DEDUP_SAMPLE            } from '../../../modules/nf-core/multiqc/main'
include { SAMTOOLS_DEPTH_MEAN   as PBQC_SAMTOOLS_DEPTH_MEAN             } from '../../../modules/local/samtools/samtools_depth_mean.nf'

workflow PROCESSED_BAM_QC {
    take:
    reference
    mq_filtered_bam
    mq_filtered_index
    rm_short_reads_bam
    rm_short_reads_bam_index
    merged_bam_lib
    merged_bam_lib_index
    dedup_lib
    dedup_lib_index
    merged_bam_sample
    merged_bam_sample_index
    dedup_sample
    dedup_sample_index
    ch_bed_file

    main:
    ch_versions                              = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Flagstat and MultiQC on processed BAM files
    ////////////////////////////////////////////////////////////////////////////////////////////////


    // mq_filtered_bam
    ch_flagstat_mq_filtered_bam              = mq_filtered_bam.join(mq_filtered_index)

    PBQC_FLAGSTAT_MQ_FILTERED_BAM ( ch_flagstat_mq_filtered_bam )
    ch_versions                              = ch_versions.mix(PBQC_FLAGSTAT_MQ_FILTERED_BAM.out.versions)

    // Run MultiQC on samtools flagstat output
    ch_multiqc_mq_filtered_bam_files         = PBQC_FLAGSTAT_MQ_FILTERED_BAM.out.flagstat.map{ meta, flagstat -> flagstat }.collect()
    ch_multiqc_config                        = params.multiqc_config       ? Channel.fromPath( params.multiqc_config,       checkIfExists: true ) : Channel.empty()
    ch_multiqc_extra_config                  = params.multiqc_extra_config ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                          = params.multiqc_logo         ? Channel.fromPath( params.multiqc_logo,         checkIfExists: true ) : Channel.empty()

    PBQC_MULTIQC_MQ_FILTERED_BAM (
        ch_multiqc_mq_filtered_bam_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(PBQC_MULTIQC_MQ_FILTERED_BAM.out.versions)

    // rm_short_reads_bam
    if (params.readlength == "auto") {

        ch_flagstat_rm_short_reads_bam           = rm_short_reads_bam.join(rm_short_reads_bam_index)

        PBQC_FLAGSTAT_RM_SHORT_READS_BAM ( ch_flagstat_rm_short_reads_bam )
        ch_versions                              = ch_versions.mix(PBQC_FLAGSTAT_RM_SHORT_READS_BAM.out.versions)

        // Run MultiQC on samtools flagstat output
        ch_multiqc_rm_short_reads_bam_files      = PBQC_FLAGSTAT_RM_SHORT_READS_BAM.out.flagstat.map{ meta, flagstat -> flagstat }.collect()
        ch_multiqc_config                        = params.multiqc_config       ? Channel.fromPath( params.multiqc_config,       checkIfExists: true ) : Channel.empty()
        ch_multiqc_extra_config                  = params.multiqc_extra_config ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true ) : Channel.empty()
        ch_multiqc_logo                          = params.multiqc_logo         ? Channel.fromPath( params.multiqc_logo,         checkIfExists: true ) : Channel.empty()

        PBQC_MULTIQC_RM_SHORT_READS_BAM (
            ch_multiqc_rm_short_reads_bam_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_extra_config.toList(),
            ch_multiqc_logo.toList()
        )
        ch_versions                              = ch_versions.mix(PBQC_MULTIQC_RM_SHORT_READS_BAM.out.versions)

    }

    // merged_bam_lib
    ch_flagstat_merged_bam_lib               = merged_bam_lib.join(merged_bam_lib_index)

    PBQC_FLAGSTAT_MERGED_BAM_LIB ( ch_flagstat_merged_bam_lib )
    ch_versions                              = ch_versions.mix(PBQC_FLAGSTAT_MERGED_BAM_LIB.out.versions)

    // Run MultiQC on samtools flagstat output
    ch_multiqc_merged_bam_lib_files          = PBQC_FLAGSTAT_MERGED_BAM_LIB.out.flagstat.map{ meta, flagstat -> flagstat }.collect()
    ch_multiqc_config                        = params.multiqc_config       ? Channel.fromPath( params.multiqc_config,       checkIfExists: true ) : Channel.empty()
    ch_multiqc_extra_config                  = params.multiqc_extra_config ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                          = params.multiqc_logo         ? Channel.fromPath( params.multiqc_logo,         checkIfExists: true ) : Channel.empty()

    PBQC_MULTIQC_MERGED_BAM_LIB (
        ch_multiqc_merged_bam_lib_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(PBQC_MULTIQC_MERGED_BAM_LIB.out.versions)

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. Preseq and Flagstat on deduplicated BAM files
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // PRESEQ
    if ( params.preseq.toBoolean() ) {
        ch_preseq_merged_bam_lib           = merged_bam_lib.join(merged_bam_lib_index)

        PBQC_PRESEQ ( ch_preseq_merged_bam_lib )
        ch_versions                              = ch_versions.mix(PBQC_PRESEQ.out.versions)

        PBQC_PLOT_PRESEQ ( PBQC_PRESEQ.out.preseq_txt )
        ch_versions                              = ch_versions.mix(PBQC_PLOT_PRESEQ.out.versions)
    }

    // flagstat dedup_lib
    ch_flagstat_dedup_lib                    = dedup_lib.join(dedup_lib_index)

    PBQC_FLAGSTAT_DEDUP_LIB ( ch_flagstat_dedup_lib )
    ch_versions                              = ch_versions.mix(PBQC_FLAGSTAT_DEDUP_LIB.out.versions)

    // Run MultiQC on samtools flagstat output
    ch_multiqc_dedup_lib_files               = PBQC_FLAGSTAT_DEDUP_LIB.out.flagstat.map{ meta, flagstat -> flagstat }.collect()

    PBQC_MULTIQC_DEDUP_LIB (
        ch_multiqc_dedup_lib_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(PBQC_MULTIQC_DEDUP_LIB.out.versions)

    // merged_bam_sample
    ch_flagstat_merged_bam_sample            = merged_bam_sample.join(merged_bam_sample_index)

    PBQC_FLAGSTAT_MERGED_BAM_SAMPLE ( ch_flagstat_merged_bam_sample )
    ch_versions                              = ch_versions.mix(PBQC_FLAGSTAT_MERGED_BAM_SAMPLE.out.versions)

    // Run MultiQC on samtools flagstat output
    ch_multiqc_merged_bam_sample_files       = PBQC_FLAGSTAT_MERGED_BAM_SAMPLE.out.flagstat.map{ meta, flagstat -> flagstat }.collect()

    PBQC_MULTIQC_MERGED_BAM_SAMPLE (
        ch_multiqc_merged_bam_sample_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(PBQC_MULTIQC_MERGED_BAM_SAMPLE.out.versions)

    // dedup_sample
    ch_dedup_sample_bam_bai                  = dedup_sample.join(dedup_sample_index)

    PBQC_FLAGSTAT_DEDUP_SAMPLE ( ch_dedup_sample_bam_bai )
    ch_versions                              = ch_versions.mix(PBQC_FLAGSTAT_DEDUP_SAMPLE.out.versions)

    // Run MultiQC on samtools flagstat output
    ch_multiqc_dedup_sample_files            = PBQC_FLAGSTAT_DEDUP_SAMPLE.out.flagstat.map{ meta, flagstat -> flagstat }.collect()

    PBQC_MULTIQC_DEDUP_SAMPLE (
        ch_multiqc_dedup_sample_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(PBQC_MULTIQC_DEDUP_SAMPLE.out.versions)

    // samtools depth mean
    PBQC_SAMTOOLS_DEPTH_MEAN (
        ch_dedup_sample_bam_bai,
        ch_bed_file
    )
    ch_versions                              = ch_versions.mix(PBQC_SAMTOOLS_DEPTH_MEAN.out.versions)

    ////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    multiqc_rm_short_reads_report            = params.readlength == "auto" ? PBQC_MULTIQC_RM_SHORT_READS_BAM.out.report.toList() : Channel.empty()      // channel: [ val(meta), path(report) ]
    multiqc_mq_filtered_report               = PBQC_MULTIQC_MQ_FILTERED_BAM.out.report.toList()                                                         // channel: [ val(meta), path(report) ]
    multiqc_merged_bam_lib_report            = PBQC_MULTIQC_MERGED_BAM_LIB.out.report.toList()                                                          // channel: [ val(meta), path(report) ]
    multiqc_dedup_lib_report                 = PBQC_MULTIQC_DEDUP_LIB.out.report.toList()                                                               // channel: [ val(meta), path(report) ]
    multiqc_merged_bam_sample_report         = PBQC_MULTIQC_MERGED_BAM_SAMPLE.out.report.toList()                                                       // channel: [ val(meta), path(report) ]
    multiqc_dedup_sample_report              = PBQC_MULTIQC_DEDUP_SAMPLE.out.report.toList()                                                            // channel: [ val(meta), path(report) ]
    dpstats                                  = PBQC_SAMTOOLS_DEPTH_MEAN.out.dpstats                                                                     // channel: [ val(meta), path(dpstats) ]
    mq_filtered_bam_flagstat                 = PBQC_FLAGSTAT_MQ_FILTERED_BAM.out.flagstat                                                               // channel: [ val(meta), path(flagstat) ]
    dedup_lib_flagstat                       = PBQC_FLAGSTAT_DEDUP_LIB.out.flagstat                                                                     // channel: [ val(meta), path(flagstat) ]
    dedup_sample_flagstat                    = PBQC_FLAGSTAT_DEDUP_SAMPLE.out.flagstat                                                                  // channel: [ val(meta), path(flagstat) ]
    preseq_txt                               = params.preseq.toBoolean() ? PBQC_PRESEQ.out.preseq_txt : Channel.empty()                                 // channel: [ val(meta), path(preseq_txt) ]
    preseq_plot                              = params.preseq.toBoolean() ? PBQC_PLOT_PRESEQ.out.preseq_plot : Channel.empty()                           // channel: [ val(meta), path(preseq_plot) ]
    versions                                 = ch_versions                                                                                                  // channel: [ versions.yml ]
}