#! /usr/bin/env nextflow

include { DEEPVARIANT_CALL       as VCD_DEEPVARIANT_CALL       } from '../../../modules/local/deepvariant/deepvariant_call.nf'
include { BCFTOOLS_STATS         as VCD_RAW_VCF_STATS          } from '../../../modules/local/bcftools/bcftools_stats.nf'
include { DEEPVARIANT_FILTER     as VCD_DEEPVARIANT_FILTER     } from '../../../modules/local/deepvariant/deepvariant_filter.nf'
include { BCFTOOLS_STATS         as VCD_FILTERED_BCF_STATS     } from '../../../modules/local/bcftools/bcftools_stats.nf'


workflow VARIANT_CALLING_DEEPVARIANT {
    take:
    bam         // channel: [ val(meta), path(bam), path(bai) ] - one entry per sample (not a joint cohort BAM list)
    reference
    fai

    main:
    ch_versions                     = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Variant calling with DeepVariant (per sample) and generating raw VCF statistics
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    VCD_DEEPVARIANT_CALL ( bam, reference, fai )
    ch_versions                     = ch_versions.mix(VCD_DEEPVARIANT_CALL.out.versions)

    VCD_RAW_VCF_STATS ( VCD_DEEPVARIANT_CALL.out.raw_vcf )
    ch_versions                     = ch_versions.mix(VCD_RAW_VCF_STATS.out.versions)


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. Filtering: depth + quality, optional SNV-only, homozygous-alt only. Generating filtered BCF statistics
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    VCD_DEEPVARIANT_FILTER ( VCD_DEEPVARIANT_CALL.out.raw_vcf )
    ch_versions                     = ch_versions.mix(VCD_DEEPVARIANT_FILTER.out.versions)

    VCD_FILTERED_BCF_STATS ( VCD_DEEPVARIANT_FILTER.out.filtered_bcf.map { meta, bcf, csi -> [ meta, bcf ] } )
    ch_versions                     = ch_versions.mix(VCD_FILTERED_BCF_STATS.out.versions)

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    raw_vcf                        = VCD_DEEPVARIANT_CALL.out.raw_vcf         // channel: [ val(meta), raw.vcf.gz ]
    raw_vcf_stats                  = VCD_RAW_VCF_STATS.out.bcf_stats          // channel: [ val(meta), raw-vcf-stats.txt ]
    filtered_bcf                   = VCD_DEEPVARIANT_FILTER.out.filtered_bcf  // channel: [ val(meta), homalt.bcf, homalt.bcf.csi ]
    filtered_bcf_stats             = VCD_FILTERED_BCF_STATS.out.bcf_stats     // channel: [ val(meta), filtered-bcf-stats.txt ]
    versions                       = ch_versions                              // channel: [ versions.yml ]
}
