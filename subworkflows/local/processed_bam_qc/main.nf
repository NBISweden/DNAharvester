#! /usr/bin/env nextflow

include { SAMTOOLS_FLAGSTAT as FLAGSTAT_MQ_FILTERED_BAM    } from '../../../modules/local/samtools/flagstat/main'
include { MULTIQC as MULTIQC_MQ_FILTERED_BAM               } from '../../../modules/nf-core/multiqc/main'
include { SAMTOOLS_FLAGSTAT as FLAGSTAT_RM_SHORT_READS_BAM } from '../../../modules/local/samtools/flagstat/main'
include { MULTIQC as MULTIQC_RM_SHORT_READS_BAM            } from '../../../modules/nf-core/multiqc/main'
include { SAMTOOLS_FLAGSTAT as FLAGSTAT_MERGED_BAM_LIB     } from '../../../modules/local/samtools/flagstat/main'
include { MULTIQC as MULTIQC_MERGED_BAM_LIB                } from '../../../modules/nf-core/multiqc/main'
include { PRESEQ                                           } from '../../../modules/local/preseq/preseq'
include { PLOT_PRESEQ                                      } from '../../../modules/local/preseq/plot_preseq'
include { SAMTOOLS_FLAGSTAT as FLAGSTAT_DEDUP_LIB          } from '../../../modules/local/samtools/flagstat/main'
include { MULTIQC as MULTIQC_DEDUP_LIB                     } from '../../../modules/nf-core/multiqc/main'
include { SAMTOOLS_FLAGSTAT as FLAGSTAT_MERGED_BAM_SAMPLE  } from '../../../modules/local/samtools/flagstat/main'
include { MULTIQC as MULTIQC_MERGED_BAM_SAMPLE             } from '../../../modules/nf-core/multiqc/main'
include { SAMTOOLS_FLAGSTAT as FLAGSTAT_DEDUP_SAMPLE       } from '../../../modules/local/samtools/flagstat/main'
include { MULTIQC as MULTIQC_DEDUP_SAMPLE                  } from '../../../modules/nf-core/multiqc/main'
include { SAMTOOLS_DEPTH_MEAN                              } from '../../../modules/local/samtools/depth_mean/main'

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
    bed

    main:
    ch_versions                              = Channel.empty()

// mq_filtered_bam
    ch_flagstat_mq_filtered_bam              = mq_filtered_bam.join(mq_filtered_index)

    FLAGSTAT_MQ_FILTERED_BAM ( ch_flagstat_mq_filtered_bam )
    ch_versions                              = ch_versions.mix(FLAGSTAT_MQ_FILTERED_BAM.out.versions)

    // Run MultiQC on samtools flagstat output
    ch_multiqc_mq_filtered_bam_files         = FLAGSTAT_MQ_FILTERED_BAM.out.flagstat.map{ meta, flagstat -> flagstat }.collect()
    ch_multiqc_config                        = params.multiqc_config       ? Channel.fromPath( params.multiqc_config,       checkIfExists: true ) : Channel.empty()
    ch_multiqc_extra_config                  = params.multiqc_extra_config ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                          = params.multiqc_logo         ? Channel.fromPath( params.multiqc_logo,         checkIfExists: true ) : Channel.empty()

    MULTIQC_MQ_FILTERED_BAM (
        ch_multiqc_mq_filtered_bam_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(MULTIQC_MQ_FILTERED_BAM.out.versions)

    // rm_short_reads_bam
    if (params.readlength == "auto") {

        ch_flagstat_rm_short_reads_bam           = rm_short_reads_bam.join(rm_short_reads_bam_index)

        FLAGSTAT_RM_SHORT_READS_BAM ( ch_flagstat_rm_short_reads_bam )
        ch_versions                              = ch_versions.mix(FLAGSTAT_RM_SHORT_READS_BAM.out.versions)

        // Run MultiQC on samtools flagstat output
        ch_multiqc_rm_short_reads_bam_files      = FLAGSTAT_RM_SHORT_READS_BAM.out.flagstat.map{ meta, flagstat -> flagstat }.collect()
        ch_multiqc_config                        = params.multiqc_config       ? Channel.fromPath( params.multiqc_config,       checkIfExists: true ) : Channel.empty()
        ch_multiqc_extra_config                  = params.multiqc_extra_config ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true ) : Channel.empty()
        ch_multiqc_logo                          = params.multiqc_logo         ? Channel.fromPath( params.multiqc_logo,         checkIfExists: true ) : Channel.empty()

        MULTIQC_RM_SHORT_READS_BAM (
            ch_multiqc_rm_short_reads_bam_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_extra_config.toList(),
            ch_multiqc_logo.toList()
        )
        ch_versions                              = ch_versions.mix(MULTIQC_RM_SHORT_READS_BAM.out.versions)

    }

    // merged_bam_lib
    ch_flagstat_merged_bam_lib               = merged_bam_lib.join(merged_bam_lib_index)

    FLAGSTAT_MERGED_BAM_LIB ( ch_flagstat_merged_bam_lib )
    ch_versions                              = ch_versions.mix(FLAGSTAT_MERGED_BAM_LIB.out.versions)

    // Run MultiQC on samtools flagstat output
    ch_multiqc_merged_bam_lib_files          = FLAGSTAT_MERGED_BAM_LIB.out.flagstat.map{ meta, flagstat -> flagstat }.collect()
    ch_multiqc_config                        = params.multiqc_config       ? Channel.fromPath( params.multiqc_config,       checkIfExists: true ) : Channel.empty()
    ch_multiqc_extra_config                  = params.multiqc_extra_config ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                          = params.multiqc_logo         ? Channel.fromPath( params.multiqc_logo,         checkIfExists: true ) : Channel.empty()

    MULTIQC_MERGED_BAM_LIB (
        ch_multiqc_merged_bam_lib_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(MULTIQC_MERGED_BAM_LIB.out.versions)

    // PRESEQ
    if ( params.preseq.toBoolean() ) {
        ch_preseq_merged_bam_lib           = merged_bam_lib.join(merged_bam_lib_index)

        PRESEQ ( ch_preseq_merged_bam_lib )
        ch_versions                              = ch_versions.mix(PRESEQ.out.versions)

        PLOT_PRESEQ ( PRESEQ.out.preseq_txt )
        ch_versions                              = ch_versions.mix(PLOT_PRESEQ.out.versions)
    }

    // flagstat dedup_lib
    ch_flagstat_dedup_lib                    = dedup_lib.join(dedup_lib_index)

    FLAGSTAT_DEDUP_LIB ( ch_flagstat_dedup_lib )
    ch_versions                              = ch_versions.mix(FLAGSTAT_DEDUP_LIB.out.versions)

    // Run MultiQC on samtools flagstat output
    ch_multiqc_dedup_lib_files               = FLAGSTAT_DEDUP_LIB.out.flagstat.map{ meta, flagstat -> flagstat }.collect()

    MULTIQC_DEDUP_LIB (
        ch_multiqc_dedup_lib_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(MULTIQC_DEDUP_LIB.out.versions)

    // merged_bam_sample
    ch_flagstat_merged_bam_sample            = merged_bam_sample.join(merged_bam_sample_index)

    FLAGSTAT_MERGED_BAM_SAMPLE ( ch_flagstat_merged_bam_sample )
    ch_versions                              = ch_versions.mix(FLAGSTAT_MERGED_BAM_SAMPLE.out.versions)

    // Run MultiQC on samtools flagstat output
    ch_multiqc_merged_bam_sample_files       = FLAGSTAT_MERGED_BAM_SAMPLE.out.flagstat.map{ meta, flagstat -> flagstat }.collect()

    MULTIQC_MERGED_BAM_SAMPLE (
        ch_multiqc_merged_bam_sample_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(MULTIQC_MERGED_BAM_SAMPLE.out.versions)

    // dedup_sample
    ch_dedup_sample_bam_bai                  = dedup_sample.join(dedup_sample_index)

    FLAGSTAT_DEDUP_SAMPLE ( ch_dedup_sample_bam_bai )
    ch_versions                              = ch_versions.mix(FLAGSTAT_DEDUP_SAMPLE.out.versions)

    // Run MultiQC on samtools flagstat output
    ch_multiqc_dedup_sample_files            = FLAGSTAT_DEDUP_SAMPLE.out.flagstat.map{ meta, flagstat -> flagstat }.collect()

    MULTIQC_DEDUP_SAMPLE (
        ch_multiqc_dedup_sample_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(MULTIQC_DEDUP_SAMPLE.out.versions)

    ch_bed_file                              = params.intervals ? Channel.fromPath(params.intervals).collect() : Channel.value([])
    // Calculate mean genome-wide depth
    SAMTOOLS_DEPTH_MEAN (
        ch_dedup_sample_bam_bai,
        ch_bed_file
    )
    ch_versions                              = ch_versions.mix(SAMTOOLS_DEPTH_MEAN.out.versions)

    emit:
    multiqc_rm_short_reads_report            = params.readlength == "auto" ? MULTIQC_RM_SHORT_READS_BAM.out.report.toList() : Channel.empty()      // channel: [ val(meta), path(report) ]
    multiqc_mq_filtered_report               = MULTIQC_MQ_FILTERED_BAM.out.report.toList()                                                         // channel: [ val(meta), path(report) ]
    multiqc_merged_bam_lib_report            = MULTIQC_MERGED_BAM_LIB.out.report.toList()                                                          // channel: [ val(meta), path(report) ]
    multiqc_dedup_lib_report                 = MULTIQC_DEDUP_LIB.out.report.toList()                                                               // channel: [ val(meta), path(report) ]
    multiqc_merged_bam_sample_report         = MULTIQC_MERGED_BAM_SAMPLE.out.report.toList()                                                       // channel: [ val(meta), path(report) ]
    multiqc_dedup_sample_report              = MULTIQC_DEDUP_SAMPLE.out.report.toList()                                                            // channel: [ val(meta), path(report) ]
    dpstats                                  = SAMTOOLS_DEPTH_MEAN.out.dpstats                                                                     // channel: [ val(meta), path(dpstats) ]
    mq_filtered_bam_flagstat                 = FLAGSTAT_MQ_FILTERED_BAM.out.flagstat                                                               // channel: [ val(meta), path(flagstat) ]
    dedup_lib_flagstat                       = FLAGSTAT_DEDUP_LIB.out.flagstat                                                                     // channel: [ val(meta), path(flagstat) ]
    preseq_txt                               = params.preseq.toBoolean() ? PRESEQ.out.preseq_txt : Channel.empty()                                 // channel: [ val(meta), path(preseq_txt) ]
    preseq_plot                              = params.preseq.toBoolean() ? PLOT_PRESEQ.out.preseq_plot : Channel.empty()                           // channel: [ val(meta), path(preseq_plot) ]
    versions                                 = ch_versions                                                                                         // channel: [ versions.yml ]
}