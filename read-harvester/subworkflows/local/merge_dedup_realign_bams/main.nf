#! /usr/bin/env nextflow

include { SAMTOOLS_FAIDX                          } from '../../../modules/local/samtools/faidx/main'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_INDEX  } from '../../../modules/local/samtools/merge/main'
include { SAMREMOVEDUP as SAMREMOVEDUP_INDEX      } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_SAMPLE } from '../../../modules/local/samtools/merge/main'
include { SAMREMOVEDUP as SAMREMOVEDUP_SAMPLE     } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_INDEX                          } from '../../../modules/nf-core/samtools/index/main'
include { PICARD_CREATESEQUENCEDICTIONARY         } from '../../../modules/local/picard/createsequencedictionary/main'
include { GATK_REALIGNERTARGETCREATOR             } from '../../../modules/local/gatk/realignertargetcreator/main'
include { GATK_INDELREALIGNER                     } from '../../../modules/local/gatk/indelrealigner/main'

workflow BAM_PROCESSING {
    take:
    reference
    bam

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_FAIDX ( reference )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    ch_bam_index_to_merge = bam.map {
        meta, bam -> [ ['id':meta.id.split("_")[0] + "_" + meta.id.split("_")[1]], bam ]
        }
        .groupTuple()

    SAMTOOLS_MERGE_INDEX ( ch_bam_index_to_merge, reference, SAMTOOLS_FAIDX.out.fai )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE_INDEX.out.versions)

    SAMREMOVEDUP_INDEX ( SAMTOOLS_MERGE_INDEX.out.bam )
    ch_versions = ch_versions.mix(SAMREMOVEDUP_INDEX.out.versions)

    ch_bam_sample_to_merge = SAMTOOLS_MERGE_INDEX.out.bam.map {
        meta, bam -> [ ['id':meta.id.split("_")[0]], bam ]
        }
        .groupTuple()

    SAMTOOLS_MERGE_SAMPLE ( ch_bam_sample_to_merge, reference, SAMTOOLS_FAIDX.out.fai )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE_SAMPLE.out.versions)

    SAMREMOVEDUP_SAMPLE ( SAMTOOLS_MERGE_SAMPLE.out.bam )
    ch_versions = ch_versions.mix(SAMREMOVEDUP_SAMPLE.out.versions)

    SAMTOOLS_INDEX ( SAMREMOVEDUP_SAMPLE.out.dedup )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    PICARD_CREATESEQUENCEDICTIONARY ( reference )
    ch_versions = ch_versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)

    ch_gatk_realignertargetcreator = SAMREMOVEDUP_SAMPLE.out.dedup.join(
        SAMTOOLS_INDEX.out.bai).groupTuple()

    GATK_REALIGNERTARGETCREATOR ( 
        ch_gatk_realignertargetcreator, 
        reference, 
        SAMTOOLS_FAIDX.out.fai, 
        PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict )
    ch_versions = ch_versions.mix(GATK_REALIGNERTARGETCREATOR.out.versions)

    ch_gatk_indelrealigner = SAMREMOVEDUP_SAMPLE.out.dedup.join(
        SAMTOOLS_INDEX.out.bai).join(
            GATK_REALIGNERTARGETCREATOR.out.intervals).groupTuple()

    GATK_INDELREALIGNER ( 
        ch_gatk_indelrealigner, 
        reference, 
        SAMTOOLS_FAIDX.out.fai, 
        PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict )
    ch_versions = ch_versions.mix(GATK_INDELREALIGNER.out.versions)

    emit:
    fai               = SAMTOOLS_FAIDX.out.fai                             // channel: path(index)
    merged_bam_index  = SAMTOOLS_MERGE_INDEX.out.bam                       // channel: [ val(meta), [ bam ] ]
    dedup_index       = SAMREMOVEDUP_INDEX.out.dedup                       // channel: [ val(meta), [ bam ] ]
    merged_bam_sample = SAMTOOLS_MERGE_SAMPLE.out.bam                      // channel: [ val(meta), [ bam ] ]
    dedup_sample      = SAMREMOVEDUP_SAMPLE.out.dedup                      // channel: [ val(meta), [ bam ] ]
    reference_dict    = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict // channel: path(reference_dict)
    realigned         = GATK_INDELREALIGNER.out.bam                        // channel: [ val(meta), [ bam ] ]
    versions          = ch_versions                                        // channel: [ versions.yml ]
}