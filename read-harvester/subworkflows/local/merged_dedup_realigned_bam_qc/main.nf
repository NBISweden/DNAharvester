#! /usr/bin/env nextflow

include { QUALIMAP_BAMQC as QUALIMAP_MERGED_BAM_INDEX  } from '../../../modules/local/qualimap/bamqc/main'
include { MULTIQC as MULTIQC_MERGED_BAM_INDEX          } from '../../../modules/nf-core/multiqc/main'
include { QUALIMAP_BAMQC as QUALIMAP_DEDUP_INDEX       } from '../../../modules/local/qualimap/bamqc/main'
include { MULTIQC as MULTIQC_DEDUP_INDEX               } from '../../../modules/nf-core/multiqc/main'
include { QUALIMAP_BAMQC as QUALIMAP_MERGED_BAM_SAMPLE } from '../../../modules/local/qualimap/bamqc/main'
include { MULTIQC as MULTIQC_MERGED_BAM_SAMPLE         } from '../../../modules/nf-core/multiqc/main'
include { QUALIMAP_BAMQC as QUALIMAP_DEDUP_SAMPLE      } from '../../../modules/local/qualimap/bamqc/main'
include { MULTIQC as MULTIQC_DEDUP_SAMPLE              } from '../../../modules/nf-core/multiqc/main'
include { QUALIMAP_BAMQC as QUALIMAP_REALIGNED         } from '../../../modules/local/qualimap/bamqc/main'
include { MULTIQC as MULTIQC_REALIGNED                 } from '../../../modules/nf-core/multiqc/main'

workflow PROCESSED_BAM_QC {
    take:
    reference
    merged_bam_index
    dedup_index
    merged_bam_sample
    dedup_sample
    realigned

    main:
    ch_versions                              = Channel.empty()

    // merged_bam_index
    QUALIMAP_MERGED_BAM_INDEX ( merged_bam_index )
    ch_versions                              = ch_versions.mix(QUALIMAP_MERGED_BAM_INDEX.out.versions)

    // Run MultiQC on QualiMap output
    ch_multiqc_merged_bam_index_files        = QUALIMAP_MERGED_BAM_INDEX.out.results.map{ meta, results -> results }.collect()
    ch_multiqc_config                        = params.multiqc_config       ? Channel.fromPath( params.multiqc_config,       checkIfExists: true ) : Channel.empty()
    ch_multiqc_extra_config                  = params.multiqc_extra_config ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                          = params.multiqc_logo         ? Channel.fromPath( params.multiqc_logo,         checkIfExists: true ) : Channel.empty()

    MULTIQC_MERGED_BAM_INDEX (
        ch_multiqc_merged_bam_index_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(MULTIQC_MERGED_BAM_INDEX.out.versions)

    // dedup_index
    QUALIMAP_DEDUP_INDEX ( dedup_index )
    ch_versions                              = ch_versions.mix(QUALIMAP_DEDUP_INDEX.out.versions)

    // Run MultiQC on QualiMap output
    ch_multiqc_dedup_index_files             = QUALIMAP_DEDUP_INDEX.out.results.map{ meta, results -> results }.collect()

    MULTIQC_DEDUP_INDEX (
        ch_multiqc_dedup_index_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(MULTIQC_DEDUP_INDEX.out.versions)

    // merged_bam_sample
    QUALIMAP_MERGED_BAM_SAMPLE ( merged_bam_sample )
    ch_versions                              = ch_versions.mix(QUALIMAP_MERGED_BAM_SAMPLE.out.versions)

    // Run MultiQC on QualiMap output
    ch_multiqc_merged_bam_sample_files       = QUALIMAP_MERGED_BAM_SAMPLE.out.results.map{ meta, results -> results }.collect()

    MULTIQC_MERGED_BAM_SAMPLE (
        ch_multiqc_merged_bam_sample_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(MULTIQC_MERGED_BAM_SAMPLE.out.versions)

    // dedup_sample
    QUALIMAP_DEDUP_SAMPLE ( dedup_sample )
    ch_versions                              = ch_versions.mix(QUALIMAP_DEDUP_SAMPLE.out.versions)

    // Run MultiQC on QualiMap output
    ch_multiqc_dedup_sample_files            = QUALIMAP_DEDUP_SAMPLE.out.results.map{ meta, results -> results }.collect()

    MULTIQC_DEDUP_SAMPLE (
        ch_multiqc_dedup_sample_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(MULTIQC_DEDUP_SAMPLE.out.versions)

    // realigned
    QUALIMAP_REALIGNED ( realigned )
    ch_versions                              = ch_versions.mix(QUALIMAP_REALIGNED.out.versions)

    // Run MultiQC on QualiMap output
    ch_multiqc_realigned_files               = QUALIMAP_REALIGNED.out.results.map{ meta, results -> results }.collect()

    MULTIQC_REALIGNED (
        ch_multiqc_realigned_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(MULTIQC_REALIGNED.out.versions)

    emit:
    multiqc_merged_bam_index_report          = MULTIQC_MERGED_BAM_INDEX.out.report.toList()       // channel: [ val(meta), path(report) ]
    multiqc_dedup_index_report               = MULTIQC_DEDUP_INDEX.out.report.toList()            // channel: [ val(meta), path(report) ]
    multiqc_merged_bam_sample_report         = MULTIQC_MERGED_BAM_SAMPLE.out.report.toList()      // channel: [ val(meta), path(report) ]
    multiqc_dedup_sample_report              = MULTIQC_DEDUP_SAMPLE.out.report.toList()           // channel: [ val(meta), path(report) ]
    multiqc_realigned_report                 = MULTIQC_REALIGNED.out.report.toList()              // channel: [ val(meta), path(report) ]
    versions                                 = ch_versions                                        // channel: [ versions.yml ]
}