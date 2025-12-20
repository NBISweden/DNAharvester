#! /usr/bin/env nextflow

// Mapping quality filter
include { SAMTOOLS_VIEW_MQ                              } from '../../../modules/local/samtools/samtools_view_mq.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_VIEW_MQ_INDEX      } from '../../../modules/nf-core/samtools/index/main'

// Read Length threshold
include { ESTIMATE_READ_LEN_CUTOFF                      } from '../../../modules/local/amber/estimate_read_len_cutoff'
include { RM_SHORT_READS                                } from '../../../modules/local/samtools/samtools_rm_short_reads.nf'
include { SAMTOOLS_INDEX as RM_SHORT_READS_INDEX        } from '../../../modules/nf-core/samtools/index/main'

// Merge BAM files per library/PCR
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_LIB          } from '../../../modules/local/samtools/samtools_merge.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_MERGE_LIB_INDEX    } from '../../../modules/nf-core/samtools/index/main'

// Remove duplicates from BAM files merged per library/PCR
include { SAMREMOVEDUP as SAMREMOVEDUP_LIB              } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_INDEX as SAMREMOVEDUP_LIB_INDEX      } from '../../../modules/nf-core/samtools/index/main'

// Running MapDamage2 (rescalling BAM if params.mapdamage2_rescale is set to true)
include { MAPDAMAGE2                                    } from '../../../modules/local/mapdamage2/main'
include { SAMTOOLS_INDEX as MAPDAMAGE2_INDEX            } from '../../../modules/nf-core/samtools/index/main'

// Removing transitions or only C to T substitutions (if params.remove_transitions and params.remove_c_to_t are set to true, respectively)
include { RM_TRANSITIONS                                } from '../../../modules/local/rm_transitions/rm_transitions.nf'
include { SAMTOOLS_INDEX as RM_TRANSITIONS_INDEX  } from '../../../modules/nf-core/samtools/index/main'

// Merge BAM files per sample
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_SAMPLE       } from '../../../modules/local/samtools/samtools_merge.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_MERGE_SAMPLE_INDEX } from '../../../modules/nf-core/samtools/index/main'

// Remove duplicates from BAM files merged per sample
include { SAMREMOVEDUP as SAMREMOVEDUP_SAMPLE           } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_INDEX as SAMREMOVEDUP_SAMPLE_INDEX   } from '../../../modules/nf-core/samtools/index/main'

// GATK Indel Realignment
include { CREATE_SEQUENCE_DICTIONARY                    } from '../../../modules/local/picard/create_sequence_dictionary.nf'
include { GATK_INDEL_REALIGNER                          } from '../../../modules/local/gatk/indel_realigner.nf'
include { SAMTOOLS_INDEX as GATK_INDEL_REALIGNER_INDEX  } from '../../../modules/nf-core/samtools/index/main'


workflow BAM_PROCESSING {
    take:
    reference
    bam
    amber_txt
    fai

    main:
    ch_versions = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. MQ filtering and read length filtering
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Filter for mapping quality (provided in custom.config)
    SAMTOOLS_VIEW_MQ ( bam )
    ch_versions = ch_versions.mix ( SAMTOOLS_VIEW_MQ.out.versions )
    SAMTOOLS_VIEW_MQ_INDEX ( SAMTOOLS_VIEW_MQ.out.bam )
    ch_versions = ch_versions.mix ( SAMTOOLS_VIEW_MQ_INDEX.out.versions )

    // Filter for minimum read length estimated from AMBER output
    if (params.readlength == "auto") {

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
            def new_meta = meta.clone()
            new_meta.id = meta.id.split("_")[0] + "_" + meta.id.split("_")[1] // remove lane info from id for merging all bams per library_id
            new_meta.remove('single_end') // remove single_end info from meta as the same library can have both single-end and paired-end data
            new_meta.remove('read_group') // remove read_group info from meta as the same library can have multiple read groups
            [ new_meta, bam ]
        }.groupTuple()

    } else {
        // Directly provide BAM channel for merging
        ch_bam_lib_to_merge = SAMTOOLS_VIEW_MQ.out.bam.map { meta, bam ->
            def new_meta = meta.clone()
            new_meta.id = meta.id.split("_")[0] + "_" + meta.id.split("_")[1] // remove lane info from id for merging all bams per library_id
            new_meta.remove('single_end') // remove single_end info from meta as the same library can have both single-end and paired-end data
            new_meta.remove('read_group') // remove read_group info from meta as the same library can have multiple read groups
            [ new_meta, bam ]
        }.groupTuple()
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. Merging BAM files per library/PCR and deduplication
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Merge BAM files per library/PCR
    SAMTOOLS_MERGE_LIB ( ch_bam_lib_to_merge, reference )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE_LIB.out.versions)
    SAMTOOLS_MERGE_LIB_INDEX ( SAMTOOLS_MERGE_LIB.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE_LIB_INDEX.out.versions)


    // Remove duplicates from BAM files merged per library/PCR
    SAMREMOVEDUP_LIB ( SAMTOOLS_MERGE_LIB.out.bam, reference )
    ch_versions = ch_versions.mix(SAMREMOVEDUP_LIB.out.versions)
    SAMREMOVEDUP_LIB_INDEX ( SAMREMOVEDUP_LIB.out.dedup )
    ch_versions = ch_versions.mix(SAMREMOVEDUP_LIB_INDEX.out.versions)
    ch_dedup_lib_bai = SAMREMOVEDUP_LIB.out.dedup.join(SAMREMOVEDUP_LIB_INDEX.out.bai)

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 3. MapDamage2 and removing transitions
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Branch reads into ancient and modern based on sample_type metadata
    ch_dedup_lib_bai.branch {
        ancient: it[0].sample_type == 'ancient'
        modern : it[0].sample_type == 'modern'
    }.set { ch_dedup_lib_bai_branched }

    // Run MapDamage2 on deduplicated BAM files merged per library/PCR (ONLY ANCIENT).
    // If params.mapdamage2_rescale is set to true, the BAM files will be rescaled.
    MAPDAMAGE2 ( ch_dedup_lib_bai_branched.ancient , reference )
    ch_versions = ch_versions.mix(MAPDAMAGE2.out.versions)

    ch_bam_processed_ancient = Channel.empty()

    // use rescaled BAM files if params.mapdamage2_rescale is set to true
    if ( params.mapdamage2_rescale.toBoolean() ) {
        MAPDAMAGE2_INDEX ( MAPDAMAGE2.out.rescaled_bam )
        ch_versions = ch_versions.mix(MAPDAMAGE2_INDEX.out.versions)
        ch_bam_processed_ancient = MAPDAMAGE2.out.rescaled_bam

    } else if ( params.remove_transitions.toBoolean() ) { // remove transitions from BAM files if params.remove_transitions is set to true
        // Remove transitions from BAM files
        RM_TRANSITIONS ( ch_dedup_lib_bai_branched.ancient )
        ch_versions = ch_versions.mix(RM_TRANSITIONS.out.versions)
        RM_TRANSITIONS_INDEX ( RM_TRANSITIONS.out.rm_trans_bam )
        ch_versions = ch_versions.mix(RM_TRANSITIONS_INDEX.out.versions)
        ch_bam_processed_ancient = RM_TRANSITIONS.out.rm_trans_bam

    } else { // if params.mapdamage2_rescale is false and params.remove_transitions is false, use deduplicated BAM files merged per library/PCR
        // Just pass through the ancient BAMs (remove bai from tuple)
        ch_bam_processed_ancient = ch_dedup_lib_bai_branched.ancient.map { meta, bam, bai -> [meta, bam] }
    }

    // Modern samples bypass MapDamage2 and RM_TRANSITIONS (remove bai from tuple)
    ch_bam_processed_modern = ch_dedup_lib_bai_branched.modern.map { meta, bam, bai -> [meta, bam] }

    // Combine processed ancient and modern BAMs and prepare for merging per sample
    ch_bam_sample_to_merge = ch_bam_processed_ancient.mix(ch_bam_processed_modern).map { meta, bam ->
        def new_meta = meta.clone()
        new_meta.id = meta.id.split("_")[0] // remove library_id for merging all bams per sample_id
        new_meta.remove('library_type') // remove library_type info from meta as the same sample can have both double-stranded and single-stranded libraries
        [ new_meta, bam ]
    }.groupTuple()

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 4. Merging BAM files per sample and deduplication
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Merge BAM files per sample
    SAMTOOLS_MERGE_SAMPLE ( ch_bam_sample_to_merge, reference )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE_SAMPLE.out.versions)
    SAMTOOLS_MERGE_SAMPLE_INDEX ( SAMTOOLS_MERGE_SAMPLE.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE_SAMPLE_INDEX.out.versions)


    // Remove duplicates from BAM files merged per sample
    SAMREMOVEDUP_SAMPLE ( SAMTOOLS_MERGE_SAMPLE.out.bam, reference )
    ch_versions = ch_versions.mix(SAMREMOVEDUP_SAMPLE.out.versions)
    SAMREMOVEDUP_SAMPLE_INDEX ( SAMREMOVEDUP_SAMPLE.out.dedup )
    ch_versions = ch_versions.mix(SAMREMOVEDUP_SAMPLE_INDEX.out.versions)
    ch_bam_deup_sample_bai = SAMREMOVEDUP_SAMPLE.out.dedup.join(SAMREMOVEDUP_SAMPLE_INDEX.out.bai)

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 5. GATK Indel Realignment
    ////////////////////////////////////////////////////////////////////////////////////////////////

    if (params.indel_realignment.toBoolean()) {
        CREATE_SEQUENCE_DICTIONARY(reference)
        ch_versions = ch_versions.mix(CREATE_SEQUENCE_DICTIONARY.out.versions)

        GATK_INDEL_REALIGNER(
            ch_bam_deup_sample_bai,
            reference,
            fai,
            CREATE_SEQUENCE_DICTIONARY.out.dict
        )
        ch_versions = ch_versions.mix(GATK_INDEL_REALIGNER.out.versions)
        GATK_INDEL_REALIGNER_INDEX ( GATK_INDEL_REALIGNER.out.realigned_bam )
        ch_versions = ch_versions.mix(GATK_INDEL_REALIGNER_INDEX.out.versions)
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    mq_filtered_bam             = SAMTOOLS_VIEW_MQ.out.bam                                                                          // channel: [ val(meta), [ bam ] ]
    mq_filtered_index           = SAMTOOLS_VIEW_MQ_INDEX.out.bai                                                                    // channel: [ val(meta), [ bai ] ]
    rm_short_reads_bam          = params.readlength == "auto" ? RM_SHORT_READS.out.bam : Channel.empty()                            // channel: [ val(meta), [ bam ] ]
    rm_short_reads_index        = params.readlength == "auto" ? RM_SHORT_READS_INDEX.out.bai : Channel.empty()                      // channel: [ val(meta), [ bai ] ]
    merged_bam_lib              = SAMTOOLS_MERGE_LIB.out.bam                                                                        // channel: [ val(meta), [ bam ] ]
    merged_bam_lib_index        = SAMTOOLS_MERGE_LIB_INDEX.out.bai                                                                  // channel: [ val(meta), [ bai ] ]
    dedup_lib                   = SAMREMOVEDUP_LIB.out.dedup                                                                        // channel: [ val(meta), [ bam ] ]
    dedup_lib_index             = SAMREMOVEDUP_LIB_INDEX.out.bai                                                                    // channel: [ val(meta), [ bai ] ]
    mapdamage2_rescaled_bam     = params.mapdamage2_rescale.toBoolean() ? MAPDAMAGE2.out.rescaled_bam : Channel.empty()             // channel: [ val(meta), [ bam ] ]
    mapdamage2_rescaled_index   = params.mapdamage2_rescale.toBoolean() ? MAPDAMAGE2_INDEX.out.bai : Channel.empty()                // channel: [ val(meta), [ bai ] ]
    rm_trans_bam                = params.remove_transitions.toBoolean() ? RM_TRANSITIONS.out.rm_trans_bam : Channel.empty()         // channel: [ val(meta), [ bam ] ]
    rm_trans_index              = params.remove_transitions.toBoolean() ? RM_TRANSITIONS_INDEX.out.bai : Channel.empty()            // channel: [ val(meta), [ bai ] ]
    merged_bam_sample           = SAMTOOLS_MERGE_SAMPLE.out.bam                                                                     // channel: [ val(meta), [ bam ] ]
    merged_bam_sample_index     = SAMTOOLS_MERGE_SAMPLE_INDEX.out.bai                                                               // channel: [ val(meta), [ bai ] ]
    dedup_sample                = SAMREMOVEDUP_SAMPLE.out.dedup                                                                     // channel: [ val(meta), [ bam ] ]
    dedup_sample_index          = SAMREMOVEDUP_SAMPLE_INDEX.out.bai                                                                 // channel: [ val(meta), [ bai ] ]
    realigned_bam               = params.indel_realignment.toBoolean() ? GATK_INDEL_REALIGNER.out.realigned_bam : Channel.empty()   // channel: [ val(meta), [ bam ] ]
    realigned_bam_index         = params.indel_realignment.toBoolean() ? GATK_INDEL_REALIGNER_INDEX.out.bai : Channel.empty()       // channel: [ val(meta), [ bai ] ]
    versions                    = ch_versions                                                                                       // channel: [ versions.yml ]
}