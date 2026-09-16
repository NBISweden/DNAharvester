#! /usr/bin/env nextflow

// Mapping quality filter
include { SAMTOOLS_VIEW_MQ              as BP_SAMTOOLS_VIEW_MQ              } from '../../../modules/local/samtools/samtools_view_mq.nf'
include { SAMTOOLS_INDEX                as BP_SAMTOOLS_VIEW_MQ_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE                as BP_SAMTOOLS_MERGE_LIB            } from '../../../modules/local/samtools/samtools_merge.nf'
include { SAMTOOLS_INDEX                as BP_SAMTOOLS_MERGE_LIB_INDEX      } from '../../../modules/nf-core/samtools/index/main'
include { SAMREMOVEDUP                  as BP_SAMREMOVEDUP_LIB              } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_INDEX                as BP_SAMREMOVEDUP_LIB_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW_SUBSAMPLE       as BP_SAMTOOLS_VIEW_SUBSAMPLE       } from '../../../modules/local/samtools/samtools_view_subsample.nf'
include { CREATE_AMBER_SAMPLESHEET      as BP_CREATE_AMBER_SAMPLESHEET      } from '../../../modules/local/amber/create_amber_samplesheet'
include { AMBER                         as BP_AMBER                         } from '../../../modules/local/amber/amber'
include { ESTIMATE_READ_LEN_CUTOFF      as BP_ESTIMATE_READ_LEN_CUTOFF      } from '../../../modules/local/amber/estimate_read_len_cutoff'
include { RM_SHORT_READS                as BP_RM_SHORT_READS                } from '../../../modules/local/samtools/samtools_rm_short_reads.nf'
include { SAMTOOLS_INDEX                as BP_RM_SHORT_READS_INDEX          } from '../../../modules/nf-core/samtools/index/main'
include { MAPDAMAGE2                    as BP_MAPDAMAGE2                    } from '../../../modules/local/mapdamage2/main'
include { SAMTOOLS_INDEX                as BP_MAPDAMAGE2_INDEX              } from '../../../modules/nf-core/samtools/index/main'
include { RM_TRANSITIONS                as BP_RM_TRANSITIONS                } from '../../../modules/local/rm_transitions/rm_transitions.nf'
include { SAMTOOLS_INDEX                as BP_RM_TRANSITIONS_INDEX          } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE                as BP_SAMTOOLS_MERGE_SAMPLE         } from '../../../modules/local/samtools/samtools_merge.nf'
include { SAMTOOLS_INDEX                as BP_SAMTOOLS_MERGE_SAMPLE_INDEX   } from '../../../modules/nf-core/samtools/index/main'
include { SAMREMOVEDUP                  as BP_SAMREMOVEDUP_SAMPLE           } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_INDEX                as BP_SAMREMOVEDUP_SAMPLE_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { CREATE_SEQUENCE_DICTIONARY    as BP_CREATE_SEQUENCE_DICTIONARY    } from '../../../modules/local/picard/create_sequence_dictionary.nf'
include { GATK_INDEL_REALIGNER          as BP_GATK_INDEL_REALIGNER          } from '../../../modules/local/gatk/indel_realigner.nf'
include { SAMTOOLS_INDEX                as BP_GATK_INDEL_REALIGNER_INDEX    } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_PROCESSING {
    take:
    reference
    bam
    fai

    main:
    ch_versions = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Merging BAM files per library/PCR; Mapping Quality Filtering; Removing duplicates
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Prepare BAM files for merging per library/PCR
    ch_bam_lib_to_merge = bam.map { meta, bam_f ->
        def new_meta = meta.clone()
        new_meta.id = meta.id.split("_")[0] + "_" + meta.id.split("_")[1] // remove lane info from id for merging all bams per library_id
        new_meta.remove('single_end') // remove single_end info from meta as the same library can have both single-end and paired-end data
        new_meta.remove('read_group') // remove read_group info from meta as the same library can have multiple read groups
        [ new_meta, bam_f ]
    }.groupTuple()

    // Merge BAM files per library/PCR
    BP_SAMTOOLS_MERGE_LIB ( ch_bam_lib_to_merge, reference )
    ch_versions = ch_versions.mix(BP_SAMTOOLS_MERGE_LIB.out.versions)
    BP_SAMTOOLS_MERGE_LIB_INDEX ( BP_SAMTOOLS_MERGE_LIB.out.bam )
    ch_versions = ch_versions.mix(BP_SAMTOOLS_MERGE_LIB_INDEX.out.versions)

    // Filter for mapping quality (provided in custom.config)
    BP_SAMTOOLS_VIEW_MQ ( BP_SAMTOOLS_MERGE_LIB.out.bam )
    ch_versions = ch_versions.mix ( BP_SAMTOOLS_VIEW_MQ.out.versions )
    BP_SAMTOOLS_VIEW_MQ_INDEX ( BP_SAMTOOLS_VIEW_MQ.out.bam )
    ch_versions = ch_versions.mix ( BP_SAMTOOLS_VIEW_MQ_INDEX.out.versions )

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. Run AMBER on MQ filtered deduplicated BAM files merged per library/PCR
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Subsample BAM files for AMBER
    BP_SAMTOOLS_VIEW_SUBSAMPLE ( BP_SAMTOOLS_VIEW_MQ.out.bam, reference )
    ch_versions = ch_versions.mix( BP_SAMTOOLS_VIEW_SUBSAMPLE.out.versions )

    // Create AMBER samplesheet
    BP_CREATE_AMBER_SAMPLESHEET ( BP_SAMTOOLS_VIEW_SUBSAMPLE.out.subsampled_bam )
    ch_versions = ch_versions.mix( BP_CREATE_AMBER_SAMPLESHEET.out.versions )

    // Run AMBER
    ch_subsampled_bam_sheet = BP_SAMTOOLS_VIEW_SUBSAMPLE.out.subsampled_bam.join( BP_CREATE_AMBER_SAMPLESHEET.out.tsv )
    BP_AMBER ( ch_subsampled_bam_sheet )
    ch_versions = ch_versions.mix( BP_AMBER.out.versions )

    // If readlength is set to "auto", estimate read length cutoff and remove short reads
    if (params.readlength == "auto") {
        // Estimate read length cutoff
        BP_ESTIMATE_READ_LEN_CUTOFF ( BP_AMBER.out.txt )
        ch_versions = ch_versions.mix ( BP_ESTIMATE_READ_LEN_CUTOFF.out.versions )

        // Remove short reads from BAM files
        ch_bam_read_len_cutoff = BP_SAMTOOLS_VIEW_MQ.out.bam.join ( BP_ESTIMATE_READ_LEN_CUTOFF.out.read_len_cutoff )
        BP_RM_SHORT_READS ( ch_bam_read_len_cutoff )
        ch_versions = ch_versions.mix ( BP_RM_SHORT_READS.out.versions )
        BP_RM_SHORT_READS_INDEX ( BP_RM_SHORT_READS.out.bam )
        ch_versions = ch_versions.mix ( BP_RM_SHORT_READS_INDEX.out.versions )

        // channel with MQ filtered and short read removed BAM files merged per library/PCR
        ch_dedup_lib_bam_input = BP_RM_SHORT_READS.out.bam
        ch_dedup_lib_bam_input_index = BP_RM_SHORT_READS_INDEX.out.bai
    }
    else {
        // othersise, just pass through the MQ filtered BAM files merged per library/PCR
        ch_dedup_lib_bam_input = BP_SAMTOOLS_VIEW_MQ.out.bam
        ch_dedup_lib_bam_input_index = BP_SAMTOOLS_VIEW_MQ_INDEX.out.bai
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 3. Removing duplicates
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Remove duplicates from BAM files merged per library/PCR
    BP_SAMREMOVEDUP_LIB ( ch_dedup_lib_bam_input, reference )
    ch_versions = ch_versions.mix(BP_SAMREMOVEDUP_LIB.out.versions)
    BP_SAMREMOVEDUP_LIB_INDEX ( BP_SAMREMOVEDUP_LIB.out.dedup )
    ch_versions = ch_versions.mix(BP_SAMREMOVEDUP_LIB_INDEX.out.versions)

    // prepare deduplicated BAM files merged per library/PCR after short read removal
    ch_dedup_lib_bai = BP_SAMREMOVEDUP_LIB.out.dedup.join(BP_SAMREMOVEDUP_LIB_INDEX.out.bai)

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 3. MapDamage2 and removing transitions on ancient samples
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Branch reads into ancient and modern based on sample_type metadata
    ch_dedup_lib_bai.branch {
        ancient: it[0].sample_type == 'ancient'
        modern : it[0].sample_type == 'modern'
    }.set { ch_dedup_lib_bai_branched }

    // Run MapDamage2 on deduplicated BAM files merged per library/PCR (ONLY ANCIENT).
    // If params.mapdamage2_rescale is set to true, the BAM files will be rescaled.
    BP_MAPDAMAGE2 ( ch_dedup_lib_bai_branched.ancient , reference )
    ch_versions = ch_versions.mix(BP_MAPDAMAGE2.out.versions)

    ch_bam_processed_ancient = Channel.empty()

    // use rescaled BAM files if params.mapdamage2_rescale is set to true
    if ( params.mapdamage2_rescale.toBoolean() ) {
        BP_MAPDAMAGE2_INDEX ( BP_MAPDAMAGE2.out.rescaled_bam )
        ch_versions = ch_versions.mix(BP_MAPDAMAGE2_INDEX.out.versions)
        ch_bam_processed_ancient = BP_MAPDAMAGE2.out.rescaled_bam

    } else if ( params.remove_transitions.toBoolean() ) { // remove transitions from BAM files if params.remove_transitions is set to true
        // Remove transitions from BAM files
        BP_RM_TRANSITIONS ( ch_dedup_lib_bai_branched.ancient )
        ch_versions = ch_versions.mix(BP_RM_TRANSITIONS.out.versions)
        BP_RM_TRANSITIONS_INDEX ( BP_RM_TRANSITIONS.out.rm_trans_bam )
        ch_versions = ch_versions.mix(BP_RM_TRANSITIONS_INDEX.out.versions)
        ch_bam_processed_ancient = BP_RM_TRANSITIONS.out.rm_trans_bam

    } else { // if params.mapdamage2_rescale is false and params.remove_transitions is false, use deduplicated BAM files merged per library/PCR
        // Just pass through the ancient BAMs (remove bai from tuple)
        ch_bam_processed_ancient = ch_dedup_lib_bai_branched.ancient.map { meta, bam_f, bai -> [meta, bam_f] }
    }

    // Modern samples bypass MapDamage2 and RM_TRANSITIONS (remove bai from tuple)
    ch_bam_processed_modern = ch_dedup_lib_bai_branched.modern.map { meta, bam_f, bai -> [meta, bam_f] }

    // Combine processed ancient and modern BAMs
    ch_bam_processed_lib = ch_bam_processed_ancient.mix(ch_bam_processed_modern)

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 4. Merging BAM files per sample and deduplication
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Combine processed ancient and modern BAMs and prepare for merging per sample
    ch_bam_sample_to_merge = ch_bam_processed_lib.map { meta, bam_f ->
        def new_meta = meta.clone()
        new_meta.id = meta.id.split("_")[0] // remove library_id for merging all bams per sample_id
        new_meta.remove('library_type') // remove library_type info from meta as the same sample can have both double-stranded and single-stranded libraries
        [ new_meta, bam_f ]
    }.groupTuple()

    // Merge BAM files per sample
    BP_SAMTOOLS_MERGE_SAMPLE ( ch_bam_sample_to_merge, reference )
    ch_versions = ch_versions.mix(BP_SAMTOOLS_MERGE_SAMPLE.out.versions)
    BP_SAMTOOLS_MERGE_SAMPLE_INDEX ( BP_SAMTOOLS_MERGE_SAMPLE.out.bam )
    ch_versions = ch_versions.mix(BP_SAMTOOLS_MERGE_SAMPLE_INDEX.out.versions)


    // Remove duplicates from BAM files merged per sample
    BP_SAMREMOVEDUP_SAMPLE ( BP_SAMTOOLS_MERGE_SAMPLE.out.bam, reference )
    ch_versions = ch_versions.mix(BP_SAMREMOVEDUP_SAMPLE.out.versions)
    BP_SAMREMOVEDUP_SAMPLE_INDEX ( BP_SAMREMOVEDUP_SAMPLE.out.dedup )
    ch_versions = ch_versions.mix(BP_SAMREMOVEDUP_SAMPLE_INDEX.out.versions)
    ch_bam_deup_sample_bai = BP_SAMREMOVEDUP_SAMPLE.out.dedup.join(BP_SAMREMOVEDUP_SAMPLE_INDEX.out.bai)

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 5. GATK Indel Realignment
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (params.indel_realignment.toBoolean()) {
        BP_CREATE_SEQUENCE_DICTIONARY(reference)
        ch_versions = ch_versions.mix(BP_CREATE_SEQUENCE_DICTIONARY.out.versions)

        BP_GATK_INDEL_REALIGNER(
            ch_bam_deup_sample_bai,
            reference,
            fai,
            BP_CREATE_SEQUENCE_DICTIONARY.out.dict.collect()
        )
        ch_versions = ch_versions.mix(BP_GATK_INDEL_REALIGNER.out.versions)
        BP_GATK_INDEL_REALIGNER_INDEX ( BP_GATK_INDEL_REALIGNER.out.realigned_bam )
        ch_versions = ch_versions.mix(BP_GATK_INDEL_REALIGNER_INDEX.out.versions)

        // GATK-realigned BAM supersedes the pre-realignment dedup BAM as "the final processed
        // sample BAM" (matches the publishDir logic below, which only publishes the
        // pre-realignment dedup BAM when indel_realignment is off)
        ch_final_dedup_sample_bam = BP_GATK_INDEL_REALIGNER.out.realigned_bam
        ch_final_dedup_sample_bai = BP_GATK_INDEL_REALIGNER_INDEX.out.bai
    } else {
        ch_final_dedup_sample_bam = BP_SAMREMOVEDUP_SAMPLE.out.dedup
        ch_final_dedup_sample_bai = BP_SAMREMOVEDUP_SAMPLE_INDEX.out.bai
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    merged_bam_lib              = BP_SAMTOOLS_MERGE_LIB.out.bam                                                                        // channel: [ val(meta), [ bam ] ]
    merged_bam_lib_index        = BP_SAMTOOLS_MERGE_LIB_INDEX.out.bai                                                                  // channel: [ val(meta), [ bai ] ]
    mq_filtered_bam             = BP_SAMTOOLS_VIEW_MQ.out.bam                                                                          // channel: [ val(meta), [ bam ] ]
    mq_filtered_index           = BP_SAMTOOLS_VIEW_MQ_INDEX.out.bai                                                                    // channel: [ val(meta), [ bai ] ]
    rm_short_reads_bam          = params.readlength == "auto" ? BP_RM_SHORT_READS.out.bam : Channel.empty()                            // channel: [ val(meta), [ bam ] ]
    rm_short_reads_index        = params.readlength == "auto" ? BP_RM_SHORT_READS_INDEX.out.bai : Channel.empty()                      // channel: [ val(meta), [ bai ] ]
    dedup_lib                   = BP_SAMREMOVEDUP_LIB.out.dedup                                                                        // channel: [ val(meta), [ bam ] ]
    dedup_lib_index             = BP_SAMREMOVEDUP_LIB_INDEX.out.bai                                                                    // channel: [ val(meta), [ bai ] ]
    mapdamage2_rescaled_bam     = params.mapdamage2_rescale.toBoolean() ? BP_MAPDAMAGE2.out.rescaled_bam : Channel.empty()             // channel: [ val(meta), [ bam ] ]
    mapdamage2_rescaled_index   = params.mapdamage2_rescale.toBoolean() ? BP_MAPDAMAGE2_INDEX.out.bai : Channel.empty()                // channel: [ val(meta), [ bai ] ]
    rm_trans_bam                = params.remove_transitions.toBoolean() ? BP_RM_TRANSITIONS.out.rm_trans_bam : Channel.empty()         // channel: [ val(meta), [ bam ] ]
    rm_trans_index              = params.remove_transitions.toBoolean() ? BP_RM_TRANSITIONS_INDEX.out.bai : Channel.empty()            // channel: [ val(meta), [ bai ] ]
    merged_bam_sample           = BP_SAMTOOLS_MERGE_SAMPLE.out.bam                                                                     // channel: [ val(meta), [ bam ] ]
    merged_bam_sample_index     = BP_SAMTOOLS_MERGE_SAMPLE_INDEX.out.bai                                                               // channel: [ val(meta), [ bai ] ]
    // "Final processed sample BAM": the GATK-realigned BAM when indel_realignment is enabled,
    // otherwise the pre-realignment dedup BAM - see ch_final_dedup_sample_bam/_bai above
    dedup_sample                = ch_final_dedup_sample_bam                                                                            // channel: [ val(meta), [ bam ] ]
    dedup_sample_index          = ch_final_dedup_sample_bai                                                                            // channel: [ val(meta), [ bai ] ]
    versions                    = ch_versions                                                                                          // channel: [ versions.yml ]
}