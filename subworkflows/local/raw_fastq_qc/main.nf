#! /usr/bin/env nextflow

include { FASTQC as FASTQC_RAW   } from '../../../modules/nf-core/fastqc/main'
include { MULTIQC as MULTIQC_RAW } from '../../../modules/nf-core/multiqc/main'

workflow RAW_FASTQ_QC {
    take:
    reads

    main:
    ch_versions = Channel.empty()

    FASTQC_RAW ( reads )
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions)

    // Run MultiQC on FastQC output
    ch_multiqc_processed_files  = FASTQC_RAW.out.zip.map{ meta, qcfile -> qcfile }.collect()
    ch_multiqc_config           = params.multiqc_config       ? Channel.fromPath( params.multiqc_config,       checkIfExists: true ) : Channel.empty()
    ch_multiqc_extra_config     = params.multiqc_extra_config ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo             = params.multiqc_logo         ? Channel.fromPath( params.multiqc_logo,         checkIfExists: true ) : Channel.empty()

    MULTIQC_RAW (
        ch_multiqc_processed_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions = ch_versions.mix(MULTIQC_RAW.out.versions)

    emit:
    fastqc_html         = FASTQC_RAW.out.html                          // channel: [ val(meta), path(html) ]
    fastqc_zip          = FASTQC_RAW.out.zip                           // channel: [ val(meta), path(zip) ]
    multiqc_report      = MULTIQC_RAW.out.report.toList()              // channel: [ val(meta), path(report) ]
    versions            = ch_versions                                  // channel: [ versions.yml ]
}
