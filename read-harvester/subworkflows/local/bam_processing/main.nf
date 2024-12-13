#! /usr/bin/env nextflow

// Mapping quality filter
include { SAMTOOLS_VIEW_MQ                              } from '../../../modules/local/samtools/view_mq/main'
include { SAMTOOLS_INDEX as SAMTOOLS_VIEW_MQ_INDEX      } from '../../../modules/nf-core/samtools/index/main'

// Read Length threshold
include { ESTIMATE_READ_LEN_CUTOFF                      } from '../../../modules/local/amber/estimate_read_len_cutoff'
include { RM_SHORT_READS                                } from '../../../modules/local/samtools/rm_short_reads/main'
include { SAMTOOLS_INDEX as RM_SHORT_READS_INDEX        } from '../../../modules/nf-core/samtools/index/main'

// Merge BAM files per library index
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_LIB          } from '../../../modules/local/samtools/merge/main'
include { SAMTOOLS_INDEX as SAMTOOLS_MERGE_LIB_INDEX    } from '../../../modules/nf-core/samtools/index/main'

// Remove duplicates from BAM files merged per library index
include { SAMREMOVEDUP as SAMREMOVEDUP_LIB              } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_INDEX as SAMREMOVEDUP_LIB_INDEX      } from '../../../modules/nf-core/samtools/index/main'

// Merge BAM files per sample
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_SAMPLE       } from '../../../modules/local/samtools/merge/main'
include { SAMTOOLS_INDEX as SAMTOOLS_MERGE_SAMPLE_INDEX } from '../../../modules/nf-core/samtools/index/main'

// Remove duplicates from BAM files merged per sample
include { SAMREMOVEDUP as SAMREMOVEDUP_SAMPLE           } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_INDEX as SAMREMOVEDUP_SAMPLE_INDEX   } from '../../../modules/nf-core/samtools/index/main'

// Realign indels
include { PICARD_CREATESEQUENCEDICTIONARY               } from '../../../modules/local/picard/createsequencedictionary/main'
include { GATK_REALIGNERTARGETCREATOR                   } from '../../../modules/local/gatk/realignertargetcreator/main'
include { GATK_INDELREALIGNER                           } from '../../../modules/local/gatk/indelrealigner/main'

workflow BAM_PROCESSING {
    take:
    reference
    fai
    bam
    bai
    amber_txt

    main:
    ch_versions = Channel.empty()

    // Create channels
    ch_reference_fai = reference.join(fai)
    ch_samtools_view_mq = bam.join(bai)

    // Filter for mapping quality (provided in custom.config)
    SAMTOOLS_VIEW_MQ ( ch_samtools_view_mq )
    ch_versions = ch_versions.mix ( SAMTOOLS_VIEW_MQ.out.versions )

    SAMTOOLS_VIEW_MQ_INDEX ( SAMTOOLS_VIEW_MQ.out.bam )
    ch_versions = ch_versions.mix ( SAMTOOLS_VIEW_MQ_INDEX.out.versions )

    // Filter for minimum read length estimated from AMBER output
    if (params.read_len_cutoff == "auto") {

        // Estimate read length cutoff
        ESTIMATE_READ_LEN_CUTOFF ( amber_txt )
        ch_versions = ch_versions.mix ( ESTIMATE_READ_LEN_CUTOFF.out.versions )

        // Remove short reads from BAM files
        ch_bam_read_len_cutoff = SAMTOOLS_VIEW_MQ.out.bam.join ( ESTIMATE_READ_LEN_CUTOFF.out.read_len_cutoff )
        RM_SHORT_READS ( ch_bam_read_len_cutoff )
        ch_versions = ch_versions.mix ( RM_SHORT_READS.out.versions )

        RM_SHORT_READS_INDEX ( RM_SHORT_READS.out.bam )
        ch_versions = ch_versions.mix ( RM_SHORT_READS_INDEX.out.versions )

        // Prepare BAM files for merging
        ch_bam_lib_to_merge = RM_SHORT_READS.out.bam.map { meta, bam ->
            [['id': meta.id.split("_")[0] + "_" + meta.id.split("_")[1]], bam]
        }.groupTuple()

    } else {
        // Directly provide BAM channel for merging
        ch_bam_lib_to_merge = SAMTOOLS_VIEW_MQ.out.bam.map { meta, bam ->
            [['id': meta.id.split("_")[0] + "_" + meta.id.split("_")[1]], bam]
        }.groupTuple()
    }

    SAMTOOLS_MERGE_LIB ( ch_bam_lib_to_merge, ch_reference_fai )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE_LIB.out.versions)

    SAMTOOLS_MERGE_LIB_INDEX ( SAMTOOLS_MERGE_LIB.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE_LIB_INDEX.out.versions)

    // Remove duplicates from BAM files merged per library index
    SAMREMOVEDUP_LIB ( SAMTOOLS_MERGE_LIB.out.bam )
    ch_versions = ch_versions.mix(SAMREMOVEDUP_LIB.out.versions)

    SAMREMOVEDUP_LIB_INDEX ( SAMREMOVEDUP_LIB.out.dedup )
    ch_versions = ch_versions.mix(SAMREMOVEDUP_LIB_INDEX.out.versions)

// Merge BAM files per sample
    ch_bam_sample_to_merge = SAMTOOLS_MERGE_LIB.out.bam.map {
        meta, bam -> [ ['id':meta.id.split("_")[0]], bam ]
        }
        .groupTuple()

    SAMTOOLS_MERGE_SAMPLE ( ch_bam_sample_to_merge, ch_reference_fai )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE_SAMPLE.out.versions)

    SAMTOOLS_MERGE_SAMPLE_INDEX ( SAMTOOLS_MERGE_SAMPLE.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE_SAMPLE_INDEX.out.versions)

// Remove duplicates from BAM files merged per sample
    SAMREMOVEDUP_SAMPLE ( SAMTOOLS_MERGE_SAMPLE.out.bam )
    ch_versions = ch_versions.mix(SAMREMOVEDUP_SAMPLE.out.versions)

    SAMREMOVEDUP_SAMPLE_INDEX ( SAMREMOVEDUP_SAMPLE.out.dedup )
    ch_versions = ch_versions.mix(SAMREMOVEDUP_SAMPLE_INDEX.out.versions)

// Realign indels
    PICARD_CREATESEQUENCEDICTIONARY ( reference )
    ch_versions = ch_versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)

    ch_gatk_realignertargetcreator = SAMREMOVEDUP_SAMPLE.out.dedup.join(
        SAMREMOVEDUP_SAMPLE_INDEX.out.bai).groupTuple()

    GATK_REALIGNERTARGETCREATOR (
        ch_gatk_realignertargetcreator,
        ch_reference_fai,
        PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict )
    ch_versions = ch_versions.mix(GATK_REALIGNERTARGETCREATOR.out.versions)

    ch_gatk_indelrealigner = SAMREMOVEDUP_SAMPLE.out.dedup.join(
        SAMREMOVEDUP_SAMPLE_INDEX.out.bai).join(
            GATK_REALIGNERTARGETCREATOR.out.intervals).groupTuple()

    GATK_INDELREALIGNER (
        ch_gatk_indelrealigner,
        ch_reference_fai,
        PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict )
    ch_versions = ch_versions.mix(GATK_INDELREALIGNER.out.versions)

    emit:
    mq_filtered_bam         = SAMTOOLS_VIEW_MQ.out.bam                                                          // channel: [ val(meta), [ bam ] ]
    mq_filtered_index       = SAMTOOLS_VIEW_MQ_INDEX.out.bai                                                    // channel: [ val(meta), [ bai ] ]
    rm_short_reads_bam      = params.read_len_cutoff == "auto" ? RM_SHORT_READS.out.bam : Channel.empty()       // channel: [ val(meta), [ bam ] ]
    rm_short_reads_index    = params.read_len_cutoff == "auto" ? RM_SHORT_READS_INDEX.out.bai : Channel.empty() // channel: [ val(meta), [ bai ] ]
    merged_bam_lib          = SAMTOOLS_MERGE_LIB.out.bam                                                        // channel: [ val(meta), [ bam ] ]
    merged_bam_lib_index    = SAMTOOLS_MERGE_LIB_INDEX.out.bai                                                  // channel: [ val(meta), [ bai ] ]
    dedup_lib               = SAMREMOVEDUP_LIB.out.dedup                                                        // channel: [ val(meta), [ bam ] ]
    dedup_lib_index         = SAMREMOVEDUP_LIB_INDEX.out.bai                                                    // channel: [ val(meta), [ bai ] ]
    merged_bam_sample       = SAMTOOLS_MERGE_SAMPLE.out.bam                                                     // channel: [ val(meta), [ bam ] ]
    merged_bam_sample_index = SAMTOOLS_MERGE_SAMPLE_INDEX.out.bai                                               // channel: [ val(meta), [ bai ] ]
    dedup_sample            = SAMREMOVEDUP_SAMPLE.out.dedup                                                     // channel: [ val(meta), [ bam ] ]
    dedup_sample_index      = SAMREMOVEDUP_SAMPLE_INDEX.out.bai                                                 // channel: [ val(meta), [ bai ] ]
    reference_dict          = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict                                // channel: path(reference_dict)
    realigned               = GATK_INDELREALIGNER.out.bam                                                       // channel: [ val(meta), [ bam, bai ] ]
    versions                = ch_versions                                                                       // channel: [ versions.yml ]
}