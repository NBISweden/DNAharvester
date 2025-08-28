#! /usr/bin/env nextflow

include { SAMTOOLS_FLAGSTAT          } from '../../../modules/local/samtools/samtools_flagstat.nf'
include { SAMTOOLS_VIEW_SUBSAMPLE    } from '../../../modules/local/samtools/samtools_view_subsample.nf'
include { CREATE_AMBER_SAMPLESHEET   } from '../../../modules/local/amber/create_amber_samplesheet'
include { AMBER                      } from '../../../modules/local/amber/amber'
include { MULTIQC as MULTIQC_BAM     } from '../../../modules/nf-core/multiqc/main'


workflow RAW_BAM_QC {
    take:
    reference
    bam             // bam file from mapping subworkflow
    bai             // bam index file

    main:
    ch_versions                              = Channel.empty()

    ch_bam_bai                               = bam.join(bai)

    SAMTOOLS_FLAGSTAT ( ch_bam_bai )
    ch_versions                              = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    // Run MultiQC samtools flagstat output
    ch_multiqc_bam_files                     = SAMTOOLS_FLAGSTAT.out.flagstat.map{ meta, flagstat -> flagstat }.collect()
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
    flagstat                                 = SAMTOOLS_FLAGSTAT.out.flagstat               // channel: [ val(meta), path(flagstat) ]
    subsampled_bam                           = SAMTOOLS_VIEW_SUBSAMPLE.out.subsampled_bam   // channel: [ val(meta), path(bam) ]
    amber_plot                               = AMBER.out.plot                               // channel: [ val(meta), path(plot) ]
    amber_txt                                = AMBER.out.txt                                // channel: [ val(meta), path(txt) ]
    multiqc_bam_report                       = MULTIQC_BAM.out.report.toList()              // channel: [ val(meta), path(report) ]
    versions                                 = ch_versions                                  // channel: [ versions.yml ]
}