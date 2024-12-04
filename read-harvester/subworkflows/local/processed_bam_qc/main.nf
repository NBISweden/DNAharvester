#! /usr/bin/env nextflow

include { SAMTOOLS_FLAGSTAT as FLAGSTAT_MERGED_BAM_LIB    } from '../../../modules/nf-core/samtools/flagstat/main'
include { MULTIQC as MULTIQC_MERGED_BAM_LIB               } from '../../../modules/nf-core/multiqc/main'
include { SAMTOOLS_FLAGSTAT as FLAGSTAT_DEDUP_LIB         } from '../../../modules/nf-core/samtools/flagstat/main'
include { MULTIQC as MULTIQC_DEDUP_LIB                    } from '../../../modules/nf-core/multiqc/main'
include { SAMTOOLS_FLAGSTAT as FLAGSTAT_MERGED_BAM_SAMPLE } from '../../../modules/nf-core/samtools/flagstat/main'
include { MULTIQC as MULTIQC_MERGED_BAM_SAMPLE            } from '../../../modules/nf-core/multiqc/main'
include { SAMTOOLS_FLAGSTAT as FLAGSTAT_DEDUP_SAMPLE      } from '../../../modules/nf-core/samtools/flagstat/main'
include { MULTIQC as MULTIQC_DEDUP_SAMPLE                 } from '../../../modules/nf-core/multiqc/main'
include { SAMTOOLS_FLAGSTAT as FLAGSTAT_REALIGNED         } from '../../../modules/nf-core/samtools/flagstat/main'
include { MULTIQC as MULTIQC_REALIGNED                    } from '../../../modules/nf-core/multiqc/main'
include { SAMTOOLS_DEPTH_MEAN                             } from '../../../modules/local/samtools/depth_mean/main'
include { PRESEQ as PRESEQ                                } from '../../../modules/local/preseq/main'

workflow PROCESSED_BAM_QC {
    take:
    reference
    merged_bam_lib
    merged_bam_lib_index
    dedup_lib
    dedup_lib_index
    merged_bam_sample
    merged_bam_sample_index
    dedup_sample
    dedup_sample_index
    realigned
    bed

    main:
    ch_versions                              = Channel.empty()

    // merged_bam_lib
    ch_flagstat_merged_bam_lib               = merged_bam_lib.join(merged_bam_lib_index)

    FLAGSTAT_MERGED_BAM_LIB ( ch_flagstat_merged_bam_lib )
    ch_versions                              = ch_versions.mix(FLAGSTAT_MERGED_BAM_LIB.out.versions)

    // RUN PRESEQ
    ch_preseq_merged_bam_lib_files           = merged_bam_lib.join(merged_bam_lib_index)
    PRESEQ ( ch_preseq_merged_bam_lib_files )
    ch_versions                              = ch_versions.mix(PRESEQ.out.versions)

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

    // dedup_lib
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
    ch_flagstat_dedup_sample                 = dedup_sample.join(dedup_sample_index)

    FLAGSTAT_DEDUP_SAMPLE ( ch_flagstat_dedup_sample )
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

    // realigned
    FLAGSTAT_REALIGNED ( realigned )
    ch_versions                              = ch_versions.mix(FLAGSTAT_REALIGNED.out.versions)

    // Run MultiQC on QualiMap output
    ch_multiqc_realigned_files               = FLAGSTAT_REALIGNED.out.flagstat.map{ meta, flagstat -> flagstat }.collect()

    MULTIQC_REALIGNED (
        ch_multiqc_realigned_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(MULTIQC_REALIGNED.out.versions)

    // Calculate mean genome-wide depth
    SAMTOOLS_DEPTH_MEAN (
        realigned,
        bed
    )
    ch_versions                              = ch_versions.mix(SAMTOOLS_DEPTH_MEAN.out.versions)

    emit:
    multiqc_merged_bam_lib_report            = MULTIQC_MERGED_BAM_LIB.out.report.toList()       // channel: [ val(meta), path(report) ]
    multiqc_dedup_lib_report                 = MULTIQC_DEDUP_LIB.out.report.toList()            // channel: [ val(meta), path(report) ]
    multiqc_merged_bam_sample_report         = MULTIQC_MERGED_BAM_SAMPLE.out.report.toList()    // channel: [ val(meta), path(report) ]
    multiqc_dedup_sample_report              = MULTIQC_DEDUP_SAMPLE.out.report.toList()         // channel: [ val(meta), path(report) ]
    multiqc_realigned_report                 = MULTIQC_REALIGNED.out.report.toList()            // channel: [ val(meta), path(report) ]
    dpstats                                  = SAMTOOLS_DEPTH_MEAN.out.dpstats                  // channel: [ val(meta), path(dpstats) ]
    dedup_lib_flagstat                       = FLAGSTAT_DEDUP_LIB.out.flagstat                  // channel: [ val(meta), path(flagstat) ]
    versions                                 = ch_versions                                      // channel: [ versions.yml ]
}