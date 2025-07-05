#! /usr/bin/env nextflow

include { BCFTOOLS_CALL                                     }   from '../../../modules/local/bcftools/bcftools_call.nf'
include { BCFTOOLS_STATS as RAW_BCF_STATS                   } from '../../../modules/local/bcftools/bcftools_stats.nf'

include { BCFTOOLS_FILTER                                   }   from '../../../modules/local/bcftools/bcftools_filter.nf'
include { BCFTOOLS_STATS as FILTERED_BCF_STATS              }   from '../../../modules/local/bcftools/bcftools_stats.nf'

include { BCFTOOLS_RM_INDELS                                }   from '../../../modules/local/bcftools/bcftools_rm_indels.nf'
include { BCFTOOLS_STATS as RM_INDELS_BCF_STATS             }   from '../../../modules/local/bcftools/bcftools_stats.nf'

include { BCFTOOLS_RM_ALLELIC_IMBALANCE                     }   from '../../../modules/local/bcftools/bcftools_rm_allelic_imbalance.nf'
include { BCFTOOLS_STATS as RM_ALLELIC_IMBALANCE_BCF_STATS  }   from '../../../modules/local/bcftools/bcftools_stats.nf'


workflow VARIANT_CALLING_BCFTOOLS {
    take:
    bam
    reference
    fai

    main:
    ch_bcf                          = Channel.empty()
    ch_versions                     = Channel.empty()


    // Variant calling with BCFtools
    BCFTOOLS_CALL ( bam, reference, fai )
    ch_bcf                          = BCFTOOLS_CALL.out.sorted_bcf
    ch_versions                     = ch_versions.mix(BCFTOOLS_CALL.out.versions)
    // BCFtools stats for raw BCF
    RAW_BCF_STATS ( ch_bcf )
    ch_versions                     = ch_versions.mix(RAW_BCF_STATS.out.versions)


    // Variant filtering
    BCFTOOLS_FILTER ( ch_bcf )
    ch_bcf                          = BCFTOOLS_FILTER.out.filtered_bcf
    ch_versions                     = ch_versions.mix(BCFTOOLS_FILTER.out.versions)
    // BCFtools stats for filtered BCF
    FILTERED_BCF_STATS ( ch_bcf )
    ch_versions                     = ch_versions.mix(FILTERED_BCF_STATS.out.versions)


    // Remove indels
    if (params.bcftools_remove_indels.toBoolean()) {
        BCFTOOLS_RM_INDELS ( ch_bcf )
        ch_bcf                      = BCFTOOLS_RM_INDELS.out.rm_indels_bcf
        ch_versions                 = ch_versions.mix(BCFTOOLS_RM_INDELS.out.versions)
        // BCFtools stats for BCF after removing indels
        RM_INDELS_BCF_STATS ( ch_bcf )
        ch_versions                 = ch_versions.mix(RM_INDELS_BCF_STATS.out.versions)
    }


    // Remove allelic imbalance
    if (params.bcftools_remove_allelic_imbalance.toBoolean()) {
        BCFTOOLS_RM_ALLELIC_IMBALANCE ( ch_bcf )
        ch_bcf                      = BCFTOOLS_RM_ALLELIC_IMBALANCE.out.rm_allelic_imbalance_bcf
        ch_versions                 = ch_versions.mix(BCFTOOLS_RM_ALLELIC_IMBALANCE.out.versions)
        // BCFtools stats for BCF after removing allelic imbalance
        RM_ALLELIC_IMBALANCE_BCF_STATS ( ch_bcf )
        ch_versions                 = ch_versions.mix(RM_ALLELIC_IMBALANCE_BCF_STATS.out.versions)
    }

    emit:
    bcftools_sorted_bcf             = BCFTOOLS_CALL.out.sorted_bcf                                                                                                          // channel: [ val(meta), sorted.bcf ]
    raw_bcf_stats                   = RAW_BCF_STATS.out.bcf_stats                                                                                                           // channel: [ val(meta), raw-bcf-stats.txt ]
    bcftools_filtered_bcf           = BCFTOOLS_FILTER.out.filtered_bcf                                                                                                      // channel: [ val(meta), filtered.bcf ]
    filtered_bcf_stats              = FILTERED_BCF_STATS.out.bcf_stats                                                                                                      // channel: [ val(meta), filtered-bcf-stats.txt ]
    bcftools_rm_indels_bcf          = params.bcftools_remove_indels.toBoolean() ? BCFTOOLS_RM_INDELS.out.rm_indels_bcf : Channel.empty()                                    // channel: [ val(meta), rm-indels.bcf ]
    rm_indels_bcf_stats             = params.bcftools_remove_indels.toBoolean() ? RM_INDELS_BCF_STATS.out.bcf_stats : Channel.empty()                                       // channel: [ val(meta), rm-indels-bcf-stats.txt ]
    bcftools_rm_allelic_imbalance   = params.bcftools_remove_allelic_imbalance.toBoolean() ? BCFTOOLS_RM_ALLELIC_IMBALANCE.out.rm_allelic_imbalance_bcf : Channel.empty()   // channel: [ val(meta), rm-allelic-imbalance.bcf ]
    rm_allelic_imbalance_bcf_stats  = params.bcftools_remove_allelic_imbalance.toBoolean() ? RM_ALLELIC_IMBALANCE_BCF_STATS.out.bcf_stats : Channel.empty()                 // channel: [ val(meta), rm-allelic-imbalance-bcf-stats.txt ]
    versions                        = ch_versions                                                                                                                           // channel: [ versions.yml ]
}