#! /usr/bin/env nextflow

include { SAMTOOLS_FLAGSTAT         as RBQC_SAMTOOLS_FLAGSTAT           } from '../../../modules/local/samtools/samtools_flagstat.nf'
include { SAMTOOLS_VIEW_SUBSAMPLE   as RBQC_SAMTOOLS_VIEW_SUBSAMPLE     } from '../../../modules/local/samtools/samtools_view_subsample.nf'
include { CREATE_AMBER_SAMPLESHEET  as RBQC_CREATE_AMBER_SAMPLESHEET    } from '../../../modules/local/amber/create_amber_samplesheet'
include { AMBER                     as RBQC_AMBER                       } from '../../../modules/local/amber/amber'
include { MULTIQC                   as RBQC_MULTIQC_BAM                 } from '../../../modules/nf-core/multiqc/main'

workflow RAW_BAM_QC {
    take:
    reference
    bam
    bai

    main:
    ch_versions                              = Channel.empty()
    ch_bam_bai                               = bam.join( bai )

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Run samtools flagstat and MultiQC
    ////////////////////////////////////////////////////////////////////////////////////////////////

    RBQC_SAMTOOLS_FLAGSTAT ( ch_bam_bai )
    ch_versions                              = ch_versions.mix( RBQC_SAMTOOLS_FLAGSTAT.out.versions )

    // Run MultiQC samtools flagstat output
    ch_multiqc_bam_files                     = RBQC_SAMTOOLS_FLAGSTAT.out.flagstat.map{ meta, flagstat -> flagstat }.collect()
    ch_multiqc_config                        = params.multiqc_config       ? Channel.fromPath( params.multiqc_config,       checkIfExists: true ) : Channel.empty()
    ch_multiqc_extra_config                  = params.multiqc_extra_config ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                          = params.multiqc_logo         ? Channel.fromPath( params.multiqc_logo,         checkIfExists: true ) : Channel.empty()

    RBQC_MULTIQC_BAM (
        ch_multiqc_bam_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix( RBQC_MULTIQC_BAM.out.versions )

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. Run AMBER
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Subsample BAM files for AMBER
    RBQC_SAMTOOLS_VIEW_SUBSAMPLE ( bam, reference )
    ch_versions                              = ch_versions.mix( RBQC_SAMTOOLS_VIEW_SUBSAMPLE.out.versions )
    // Create AMBER samplesheet
    RBQC_CREATE_AMBER_SAMPLESHEET ( RBQC_SAMTOOLS_VIEW_SUBSAMPLE.out.subsampled_bam )
    ch_versions                              = ch_versions.mix( RBQC_CREATE_AMBER_SAMPLESHEET.out.versions )
    // Run AMBER
    RBQC_AMBER ( RBQC_SAMTOOLS_VIEW_SUBSAMPLE.out.subsampled_bam.join( RBQC_CREATE_AMBER_SAMPLESHEET.out.tsv ))
    ch_versions                              = ch_versions.mix( RBQC_AMBER.out.versions )

    ////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    flagstat                                 = RBQC_SAMTOOLS_FLAGSTAT.out.flagstat               // channel: [ val(meta), path(flagstat) ]
    subsampled_bam                           = RBQC_SAMTOOLS_VIEW_SUBSAMPLE.out.subsampled_bam   // channel: [ val(meta), path(bam) ]
    amber_plot                               = RBQC_AMBER.out.plot                               // channel: [ val(meta), path(plot) ]
    amber_txt                                = RBQC_AMBER.out.txt                                // channel: [ val(meta), path(txt) ]
    multiqc_bam_report                       = RBQC_MULTIQC_BAM.out.report.toList()              // channel: [ val(meta), path(report) ]
    versions                                 = ch_versions                                       // channel: [ versions.yml ]
}