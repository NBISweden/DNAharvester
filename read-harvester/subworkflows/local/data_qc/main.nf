#! /usr/bin/env nextflow

include { FASTQC as FASTQC_RAW   } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_MERGE } from '../../../modules/nf-core/fastqc/main'
include { MAPDAMAGE2             } from '../../../modules/nf-core/mapdamage2/main'
include { MULTIQC                } from '../../../modules/nf-core/multiqc/main'

workflow DATA_QC {
    take:
    reference
    raw_reads     // paired-end reads or single-end reads
    trimmed_reads // merged paired-end reads or trimmed single-end reads
    bam

    main:
    ch_versions = Channel.empty()

    FASTQC_RAW ( raw_reads )
    FASTQC_TRIM ( trimmed_reads )
    ch_versions = ch_versions.mix(FASTQC_TRIM.out.versions)

    //MAPDAMAGE2 ( BWA_SAMSE.out.bam, reference )
    //ch_versions = ch_versions.mix(MAPDAMAGE2.out.versions)

    //MULTIQC (  )
    //ch_versions = ch_versions.mix(MULTIQC.out.versions)

    emit:
    fastqc_raw_html  = FASTQC_TRIM.out.html                    // channel: [ val(meta), path(html) ]
    fastqc_raw_zip   = FASTQC_TRIM.out.zip                     // channel: [ val(meta), path(zip) ]
    fastqc_trim_html = FASTQC_TRIM.out.html                    // channel: [ val(meta), path(html) ]
    fastqc_trim_zip  = FASTQC_TRIM.out.zip                     // channel: [ val(meta), path(zip) ]
    versions         = ch_versions                             // channel: [ versions.yml ]
}