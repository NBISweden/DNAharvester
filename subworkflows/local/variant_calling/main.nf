#! /usr/bin/env nextflow

include { ANGSD_VARIANT_CALLING     }       from '../../../modules/local/angsd/variant_calling/main'
include { BCFTOOLS_VARIANT_CALLING  }       from '../../../modules/local/bcftools/variant_calling/main'

workflow VARIANT_CALLING {
    take:
    bam
    reference
    fai

    main:
    ch_versions             = Channel.empty()

    ANGSD_VARIANT_CALLING ( bam, reference, fai )
    ch_versions             = ch_versions.mix(ANGSD_VARIANT_CALLING.out.versions)

    BCFTOOLS_VARIANT_CALLING ( bam, reference, fai )
    ch_versions             = ch_versions.mix(BCFTOOLS_VARIANT_CALLING.out.versions)


    emit:
    angsd_log               = ANGSD_VARIANT_CALLING.out.angsd_log               // channel: [ val(meta), log ]
    angsd_bamlist           = ANGSD_VARIANT_CALLING.out.angsd_bamlist           // channel: [ val(meta), bamlist ]
    angsd_geno              = ANGSD_VARIANT_CALLING.out.angsd_geno              // channel: [ val(meta), geno ]
    angsd_mafs              = ANGSD_VARIANT_CALLING.out.angsd_mafs              // channel: [ val(meta), mafs ]
    angsd_beagle            = ANGSD_VARIANT_CALLING.out.angsd_beagle            // channel: [ val(meta), beagle ]
    angsd_bcf               = ANGSD_VARIANT_CALLING.out.angsd_bcf               // channel: [ val(meta), bcf ]
    bcftools_sorted_bcf     = BCFTOOLS_VARIANT_CALLING.out.bcftools_sorted_bcf  // channel: [ val(meta), sorted.bcf ]
    bcftools_bcf_stats      = BCFTOOLS_VARIANT_CALLING.out.bcftools_bcf_stats   // channel: [ val(meta), bcf.stats ]
    versions                = ch_versions                                       // channel: [ versions.yml ]
}