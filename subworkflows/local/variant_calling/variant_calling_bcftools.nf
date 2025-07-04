#! /usr/bin/env nextflow

include { BCFTOOLS_CALL                 }   from '../../../modules/local/bcftools/bcftools_call.nf'
include { BCFTOOLS_FILTER               }   from '../../../modules/local/bcftools/bcftools_filter.nf'
include { BCFTOOLS_RM_INDELS            }   from '../../../modules/local/bcftools/bcftools_rm_indels.nf'
include { BCFTOOLS_RM_ALLELIC_IMBALANCE }   from '../../../modules/local/bcftools/bcftools_rm_allelic_imbalance.nf'
include { BCFTOOLS_STATS                }   from '../../../modules/local/bcftools/bcftools_stats.nf'


workflow VARIANT_CALLING_BCFTOOLS {
    take:
    bam
    reference
    fai

    main:
    ch_versions                     = Channel.empty()

    // Variant calling with BCFtools
    BCFTOOLS_CALL ( bam, reference, fai )
    ch_versions                     = ch_versions.mix(BCFTOOLS_CALL.out.versions)
    ch_bcf                          = BCFTOOLS_CALL.out.sorted_bcf

    // Variant filtering
    BCFTOOLS_FILTER ( ch_bcf )
    ch_bcf                          = BCFTOOLS_FILTER.out.filtered_bcf
    ch_versions                     = ch_versions.mix(BCFTOOLS_FILTER.out.versions)

    // Remove indels
    if (params.remove_indels.toBoolean()) {
        BCFTOOLS_RM_INDELS ( ch_bcf )
        ch_bcf                      = BCFTOOLS_RM_INDELS.out.rm_indels_bcf
        ch_versions                 = ch_versions.mix(BCFTOOLS_RM_INDELS.out.versions)
    }

    // Remove allelic imbalance
    if (params.remove_allelic_imbalance.toBoolean()) {
        BCFTOOLS_RM_ALLELIC_IMBALANCE ( ch_bcf )
        ch_bcf                      = BCFTOOLS_RM_ALLELIC_IMBALANCE.out.rm_allelic_imbalance_bcf
        ch_versions                 = ch_versions.mix(BCFTOOLS_RM_ALLELIC_IMBALANCE.out.versions)
    }

    // Generate BCFtools stats
    BCFTOOLS_STATS ( ch_bcf )
    ch_versions         = ch_versions.mix(BCFTOOLS_STATS.out.versions)

    emit:
    bcftools_sorted_bcf             = BCFTOOLS_CALL.out.sorted_bcf                                                                                      // channel: [ val(meta), sorted.bcf ]
    bcftools_filtered_bcf           = BCFTOOLS_FILTER.out.filtered_bcf                                                                                  // channel: [ val(meta), filtered.bcf ]
    bcftools_rm_indels_bcf          = params.remove_indels.toBoolean() ? BCFTOOLS_RM_INDELS.out.rm_indels_bcf : Channel.empty()                                     // channel: [ val(meta), rm-indels.bcf ]
    bcftools_rm_allelic_imbalance   = params.remove_allelic_imbalance.toBoolean() ? BCFTOOLS_RM_ALLELIC_IMBALANCE.out.rm_allelic_imbalance_bcf : Channel.empty()    // channel: [ val(meta), rm-allelic-imbalance.bcf ]
    versions                        = ch_versions                                                                                                       // channel: [ versions.yml ]
}