#! /usr/bin/env nextflow

include { FASTQC as FASTQC_RAW       } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_PROCESSED } from '../../../modules/nf-core/fastqc/main'
include { MAPDAMAGE2                 } from '../../../modules/nf-core/mapdamage2/main'
include { MULTIQC                    } from '../../../modules/nf-core/multiqc/main'

workflow DATA_QC {
    take:
    reference
    raw_reads       // paired-end reads or single-end reads
    processed_reads // merged paired-end reads or trimmed single-end reads
    bam

    main:
    ch_versions = Channel.empty()

    FASTQC_RAW ( raw_reads )
    FASTQC_PROCESSED ( processed_reads )
    ch_versions = ch_versions.mix(FASTQC_PROCESSED.out.versions)

    //MAPDAMAGE2 ( bam, reference )
    //ch_versions = ch_versions.mix(MAPDAMAGE2.out.versions)

    //MULTIQC (  )
    //ch_versions = ch_versions.mix(MULTIQC.out.versions)

    emit:
    fastqc_raw_html       = FASTQC_RAW.out.html                        // channel: [ val(meta), path(html) ]
    fastqc_raw_zip        = FASTQC_RAW.out.zip                         // channel: [ val(meta), path(zip) ]
    fastqc_processed_html = FASTQC_PROCESSED.out.html                  // channel: [ val(meta), path(html) ]
    fastqc_processed_zip  = FASTQC_PROCESSED.out.zip                   // channel: [ val(meta), path(zip) ]
    versions              = ch_versions                                // channel: [ versions.yml ]
}