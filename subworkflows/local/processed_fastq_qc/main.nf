#! /usr/bin/env nextflow

include { FASTQC as FASTQC_PROCESSED } from '../../../modules/nf-core/fastqc/main'
include { MULTIQC as MULTIQC_FASTQ   } from '../../../modules/nf-core/multiqc/main'


workflow PROCESSED_FASTQ_QC {
    take:
    processed_reads // merged and filtered paired-end reads
    fastp_json      // read statistic files from FastP

    main:
    ch_versions                              = Channel.empty()

    FASTQC_PROCESSED ( processed_reads )
    ch_versions                              = ch_versions.mix(FASTQC_PROCESSED.out.versions)

    // Run MultiQC on FastQC output
    ch_multiqc_fastq_files                   = FASTQC_PROCESSED.out.zip.map{ meta, qcfile -> qcfile }.mix(
                                                fastp_json.map{ meta, fastp_json -> fastp_json }).collect()
    ch_multiqc_config                        = params.multiqc_config       ? Channel.fromPath( params.multiqc_config,       checkIfExists: true ) : Channel.empty()
    ch_multiqc_extra_config                  = params.multiqc_extra_config ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                          = params.multiqc_logo         ? Channel.fromPath( params.multiqc_logo,         checkIfExists: true ) : Channel.empty()

    MULTIQC_FASTQ (
        ch_multiqc_fastq_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(MULTIQC_FASTQ.out.versions)


    emit:
    fastqc_html                              = FASTQC_PROCESSED.out.html                    // channel: [ val(meta), path(html) ]
    fastqc_zip                               = FASTQC_PROCESSED.out.zip                     // channel: [ val(meta), path(zip) ]
    multiqc_fastq_report                     = MULTIQC_FASTQ.out.report.toList()            // channel: [ val(meta), path(report) ]
    versions                                 = ch_versions                                  // channel: [ versions.yml ]
}