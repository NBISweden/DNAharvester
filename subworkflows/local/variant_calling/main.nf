#! /usr/bin/env nextflow

include { ANGSD_VARIANT_CALLING     }       from '../../../modules/local/angsd/variant_calling/main'
include { BCFTOOLS_CALL  }       from '../../../modules/local/bcftools/bcftools_call.nf'
include { BCFTOOLS_FILTER }      from '../../../modules/local/bcftools/bcftools_filter.nf'
include { BCFTOOLS_RM_INDELS }     from '../../../modules/local/bcftools/bcftools_rm_indels.nf'

include { BCFTOOLS_STATS } from '../../../modules/local/bcftools/bcftools_stats.nf'


workflow VARIANT_CALLING {
    take:
    bam
    reference
    fai

    main:
    ch_versions             = Channel.empty()

    // Variant calling with ANGSD
    ANGSD_VARIANT_CALLING ( bam, reference, fai )
    ch_versions             = ch_versions.mix(ANGSD_VARIANT_CALLING.out.versions)

    // Variant calling with BCFtools
    BCFTOOLS_CALL ( bam, reference, fai )
    ch_versions             = ch_versions.mix(BCFTOOLS_CALL.out.versions)

    ch_bcf = BCFTOOLS_CALL.out.sorted_bcf

    // Filter Variants for bcf files produced by BCFtools
    BCFTOOLS_FILTER ( ch_bcf )
    ch_bcf = BCFTOOLS_FILTER.out.filtered_bcf
    ch_versions             = ch_versions.mix(BCFTOOLS_FILTER.out.versions)

    if (params.remove_indels) {
        // Remove indels from BCFtools filtered variants
        BCFTOOLS_RM_INDELS ( ch_bcf )
        ch_bcf             = BCFTOOLS_RM_INDELS.out.rm_indels_bcf
        ch_versions         = ch_versions.mix(BCFTOOLS_RM_INDELS.out.versions)

    }
    if (params.remove_allelic_imbalance) {
        // Remove allelic imbalance from BCFtools filtered variants
        BCFTOOLS_RM_ALLELIC_IMBALANCE ( ch_bcf )
        ch_bcf             = BCFTOOLS_RM_ALLELIC_IMBALANCE.out.rmindels_bcf
        ch_versions         = ch_versions.mix(BCFTOOLS_RM_ALLELIC_IMBALANCE.out.versions)
    }




    emit:
    angsd_log               = ANGSD_VARIANT_CALLING.out.angsd_log               // channel: [ val(meta), log ]
    angsd_bamlist           = ANGSD_VARIANT_CALLING.out.angsd_bamlist           // channel: [ val(meta), bamlist ]
    angsd_geno              = ANGSD_VARIANT_CALLING.out.angsd_geno              // channel: [ val(meta), geno ]
    angsd_mafs              = ANGSD_VARIANT_CALLING.out.angsd_mafs              // channel: [ val(meta), mafs ]
    angsd_beagle            = ANGSD_VARIANT_CALLING.out.angsd_beagle            // channel: [ val(meta), beagle ]
    angsd_bcf               = ANGSD_VARIANT_CALLING.out.angsd_bcf               // channel: [ val(meta), bcf ]
    bcftools_sorted_bcf     = BCFTOOLS_CALL.out.sorted_bcf                      // channel: [ val(meta), sorted.bcf ]
    bcftools_filtered_bcf   = BCFTOOLS_FILTER.out.filtered_bcf                  // channel: [ val(meta), filtered.bcf ]
    versions                = ch_versions                                       // channel: [ versions.yml ]
}