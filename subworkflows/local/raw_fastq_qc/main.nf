#! /usr/bin/env nextflow

include { FASTQC    as RFQC_FASTQC      } from '../../../modules/nf-core/fastqc/main'
include { MULTIQC   as RFQC_MULTIQC     } from '../../../modules/nf-core/multiqc/main'

workflow RAW_FASTQ_QC {
    take:
    reads

    main:
    ch_versions = Channel.empty()

    RFQC_FASTQC ( reads )
    ch_versions = ch_versions.mix(RFQC_FASTQC.out.versions)

    // Run MultiQC on FastQC output
    ch_multiqc_processed_files  = RFQC_FASTQC.out.zip.map{ meta, qcfile -> qcfile }.collect()
    ch_multiqc_config           = params.multiqc_config       ? Channel.fromPath( params.multiqc_config,       checkIfExists: true ) : Channel.empty()
    ch_multiqc_extra_config     = params.multiqc_extra_config ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo             = params.multiqc_logo         ? Channel.fromPath( params.multiqc_logo,         checkIfExists: true ) : Channel.empty()

    RFQC_MULTIQC (
        ch_multiqc_processed_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions = ch_versions.mix( RFQC_MULTIQC.out.versions )

    emit:
    fastqc_html         = RFQC_FASTQC.out.html                  // channel: [ val(meta), path(html) ]
    fastqc_zip          = RFQC_FASTQC.out.zip                   // channel: [ val(meta), path(zip) ]
    multiqc_report      = RFQC_MULTIQC.out.report.toList()      // channel: [ val(meta), path(report) ]
    versions            = ch_versions                           // channel: [ versions.yml ]
}
