#! /usr/bin/env nextflow

include { FASTQC    as PFQC_FASTQC       } from '../../../modules/nf-core/fastqc/main'
include { MULTIQC   as PFQC_MULTIQC      } from '../../../modules/nf-core/multiqc/main'

workflow PROCESSED_FASTQ_QC {
    take:
    processed_reads
    fastp_json

    main:
    ch_versions                 = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Run FastQC and MultiQC on processed FASTQ files
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    PFQC_FASTQC ( processed_reads )
    ch_versions                 = ch_versions.mix(PFQC_FASTQC.out.versions)

    // Run MultiQC on FastQC output
    ch_multiqc_processed_files  = PFQC_FASTQC.out.zip.map{ meta, qcfile -> qcfile }.mix(
                                    fastp_json.map{ meta, json -> json }).collect()
    ch_multiqc_config           = params.multiqc_config       ? Channel.fromPath( params.multiqc_config,       checkIfExists: true ) : Channel.empty()
    ch_multiqc_extra_config     = params.multiqc_extra_config ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo             = params.multiqc_logo         ? Channel.fromPath( params.multiqc_logo,         checkIfExists: true ) : Channel.empty()

    PFQC_MULTIQC (
        ch_multiqc_processed_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                 = ch_versions.mix(PFQC_MULTIQC.out.versions)

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    fastqc_html                 = PFQC_FASTQC.out.html                    // channel: [ val(meta), path(html) ]
    fastqc_zip                  = PFQC_FASTQC.out.zip                     // channel: [ val(meta), path(zip) ]
    multiqc_processed_report    = PFQC_MULTIQC.out.report.toList()        // channel: [ val(meta), path(report) ]
    versions                    = ch_versions                             // channel: [ versions.yml ]
}