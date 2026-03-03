#! /usr/bin/env nextflow

include { BCFTOOLS_CALL                     as VCB_BCFTOOLS_CALL                       } from '../../../modules/local/bcftools/bcftools_call.nf'
include { BCFTOOLS_STATS                    as VCB_RAW_BCF_STATS                       } from '../../../modules/local/bcftools/bcftools_stats.nf'
include { BCFTOOLS_FILTER                   as VCB_BCFTOOLS_FILTER                     } from '../../../modules/local/bcftools/bcftools_filter.nf'
include { BCFTOOLS_STATS                    as VCB_FILTERED_BCF_STATS                  } from '../../../modules/local/bcftools/bcftools_stats.nf'
include { BCFTOOLS_RM_INDELS                as VCB_BCFTOOLS_RM_INDELS                  } from '../../../modules/local/bcftools/bcftools_rm_indels.nf'
include { BCFTOOLS_STATS                    as VCB_RM_INDELS_BCF_STATS                 } from '../../../modules/local/bcftools/bcftools_stats.nf'
include { BCFTOOLS_RM_ALLELIC_IMBALANCE     as VCB_BCFTOOLS_RM_ALLELIC_IMBALANCE       } from '../../../modules/local/bcftools/bcftools_rm_allelic_imbalance.nf'
include { BCFTOOLS_STATS                    as VCB_RM_ALLELIC_IMBALANCE_BCF_STATS      } from '../../../modules/local/bcftools/bcftools_stats.nf'


workflow VARIANT_CALLING_BCFTOOLS {
    take:
    bam
    reference
    fai
    regionsfile

    main:
    ch_versions                     = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Variant calling with BCFtools and generating raw BCF statistics
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    VCB_BCFTOOLS_CALL ( bam, reference, fai, regionsfile )
    ch_bcf                          = VCB_BCFTOOLS_CALL.out.sorted_bcf
    ch_versions                     = ch_versions.mix(VCB_BCFTOOLS_CALL.out.versions)

    VCB_RAW_BCF_STATS ( ch_bcf )
    ch_versions                     = ch_versions.mix(VCB_RAW_BCF_STATS.out.versions)


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. Variant filtering and generating filtered BCF statistics
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    VCB_BCFTOOLS_FILTER ( ch_bcf )
    ch_bcf                          = VCB_BCFTOOLS_FILTER.out.filtered_bcf
    ch_versions                     = ch_versions.mix(VCB_BCFTOOLS_FILTER.out.versions)

    VCB_FILTERED_BCF_STATS ( ch_bcf )
    ch_versions                     = ch_versions.mix(VCB_FILTERED_BCF_STATS.out.versions)


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 3. Remove indels (optional) and generate corresponding statistics
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (params.bcftools_remove_indels.toBoolean()) {
        VCB_BCFTOOLS_RM_INDELS ( ch_bcf )
        ch_bcf                      = VCB_BCFTOOLS_RM_INDELS.out.rm_indels_bcf
        ch_versions                 = ch_versions.mix(VCB_BCFTOOLS_RM_INDELS.out.versions)

        VCB_RM_INDELS_BCF_STATS ( ch_bcf )
        ch_versions                 = ch_versions.mix(VCB_RM_INDELS_BCF_STATS.out.versions)
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 4. Remove allelic imbalance (optional) and generate corresponding statistics
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (params.bcftools_remove_allelic_imbalance.toBoolean()) {
        VCB_BCFTOOLS_RM_ALLELIC_IMBALANCE ( ch_bcf )
        ch_bcf                      = VCB_BCFTOOLS_RM_ALLELIC_IMBALANCE.out.rm_allelic_imbalance_bcf
        ch_versions                 = ch_versions.mix(VCB_BCFTOOLS_RM_ALLELIC_IMBALANCE.out.versions)

        VCB_RM_ALLELIC_IMBALANCE_BCF_STATS ( ch_bcf )
        ch_versions                 = ch_versions.mix(VCB_RM_ALLELIC_IMBALANCE_BCF_STATS.out.versions)
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    bcftools_sorted_bcf             = VCB_BCFTOOLS_CALL.out.sorted_bcf                                                                                                          // channel: [ val(meta), sorted.bcf ]
    raw_bcf_stats                   = VCB_RAW_BCF_STATS.out.bcf_stats                                                                                                           // channel: [ val(meta), raw-bcf-stats.txt ]
    bcftools_filtered_bcf           = VCB_BCFTOOLS_FILTER.out.filtered_bcf                                                                                                      // channel: [ val(meta), filtered.bcf ]
    filtered_bcf_stats              = VCB_FILTERED_BCF_STATS.out.bcf_stats                                                                                                      // channel: [ val(meta), filtered-bcf-stats.txt ]
    bcftools_rm_indels_bcf          = params.bcftools_remove_indels.toBoolean() ? VCB_BCFTOOLS_RM_INDELS.out.rm_indels_bcf : Channel.empty()                                    // channel: [ val(meta), rm-indels.bcf ]
    rm_indels_bcf_stats             = params.bcftools_remove_indels.toBoolean() ? VCB_RM_INDELS_BCF_STATS.out.bcf_stats : Channel.empty()                                       // channel: [ val(meta), rm-indels-bcf-stats.txt ]
    bcftools_rm_allelic_imbalance   = params.bcftools_remove_allelic_imbalance.toBoolean() ? VCB_BCFTOOLS_RM_ALLELIC_IMBALANCE.out.rm_allelic_imbalance_bcf : Channel.empty()   // channel: [ val(meta), rm-allelic-imbalance.bcf ]
    rm_allelic_imbalance_bcf_stats  = params.bcftools_remove_allelic_imbalance.toBoolean() ? VCB_RM_ALLELIC_IMBALANCE_BCF_STATS.out.bcf_stats : Channel.empty()                 // channel: [ val(meta), rm-allelic-imbalance-bcf-stats.txt ]
    versions                        = ch_versions                                                                                                                               // channel: [ versions.yml ]
}