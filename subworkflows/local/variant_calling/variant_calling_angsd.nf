#! /usr/bin/env nextflow

include { ANGSD_VARIANT_CALLING     as VCA_ANGSD_VARIANT_CALLING    } from '../../../modules/local/angsd/variant_calling/main'
include { BCFTOOLS_STATS            as VCA_ANGSD_RAW_BCF_STATS      } from '../../../modules/local/bcftools/bcftools_stats.nf'


workflow VARIANT_CALLING_ANGSD {
    take:
    bam
    reference
    fai

    main:
    ch_versions                     = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Variant calling with ANGSD and generating raw BCF statistics
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    VCA_ANGSD_VARIANT_CALLING ( bam, reference, fai )
    ch_bcf                          = VCA_ANGSD_VARIANT_CALLING.out.angsd_bcf
    ch_versions                     = ch_versions.mix(VCA_ANGSD_VARIANT_CALLING.out.versions)

    VCA_ANGSD_RAW_BCF_STATS ( ch_bcf )
    ch_versions                     = ch_versions.mix(VCA_ANGSD_RAW_BCF_STATS.out.versions)

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    emit:
    angsd_log               = VCA_ANGSD_VARIANT_CALLING.out.angsd_log           // channel: [ val(meta), log ]
    angsd_bamlist           = VCA_ANGSD_VARIANT_CALLING.out.angsd_bamlist       // channel: [ val(meta), bamlist ]
    angsd_geno              = VCA_ANGSD_VARIANT_CALLING.out.angsd_geno          // channel: [ val(meta), geno ]
    angsd_mafs              = VCA_ANGSD_VARIANT_CALLING.out.angsd_mafs          // channel: [ val(meta), mafs ]
    angsd_beagle            = VCA_ANGSD_VARIANT_CALLING.out.angsd_beagle        // channel: [ val(meta), beagle ]
    angsd_bcf               = VCA_ANGSD_VARIANT_CALLING.out.angsd_bcf           // channel: [ val(meta), bcf ]
    angsd_raw_bcf_stats     = VCA_ANGSD_RAW_BCF_STATS.out.bcf_stats             // channel: [ val(meta), angsd_raw-bcf-stats.txt ]
    versions                = ch_versions                                       // channel: [ versions.yml ]
}