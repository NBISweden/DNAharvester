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
include { GATK_INDEL_REALIGNER                          } from '../../../modules/local/gatk/indel_realigner.nf'


workflow BAM_PROCESSING {
    take:
    reference
    bam
    amber_txt

    main:
    ch_versions = Channel.empty()


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
            // update only the 'id' field in meta, keep all other fields
            [['id': meta.id.split("_")[0] + "_" + meta.id.split("_")[1], 'library_type': meta.library_type, 'single_end': meta.single_end], bam]
        }.groupTuple()

    } else {
        // Directly provide BAM channel for merging
        ch_bam_lib_to_merge = SAMTOOLS_VIEW_MQ.out.bam.map { meta, bam ->
            // update only the 'id' field in meta, keep all other fields
            [['id': meta.id.split("_")[0] + "_" + meta.id.split("_")[1], 'library_type': meta.library_type, 'single_end': meta.single_end], bam]
        }.groupTuple()
    }


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


    // Run MapDamage2 on deduplicated BAM files merged per library/PCR. If params.mapdamage2_rescale is set to true, the BAM files will be rescaled.
    MAPDAMAGE2 ( ch_dedup_lib_bai , reference )
    ch_versions = ch_versions.mix(MAPDAMAGE2.out.versions)


    // use rescaled BAM files if params.mapdamage2_rescale is set to true
    if ( params.mapdamage2_rescale.toBoolean() ) {
        MAPDAMAGE2_INDEX ( MAPDAMAGE2.out.rescaled_bam )
        ch_versions = ch_versions.mix(MAPDAMAGE2_INDEX.out.versions)

        // Prepare BAM files for merging per sample
        ch_bam_sample_to_merge = MAPDAMAGE2.out.rescaled_bam.map { meta, bam ->
            // update only the 'id' field in meta, keep all other fields
            [['id': meta.id.split("_")[0]], bam]
        }.groupTuple()
    } else if ( params.remove_transitions.toBoolean() ) { // remove transitions from BAM files if params.remove_transitions is set to true
        // Remove transitions from BAM files
        RM_TRANSITIONS ( ch_dedup_lib_bai )
        ch_versions = ch_versions.mix(RM_TRANSITIONS.out.versions)
        RM_TRANSITIONS_INDEX ( RM_TRANSITIONS.out.rm_trans_bam,  )
        ch_versions = ch_versions.mix(RM_TRANSITIONS_INDEX.out.versions)

        // Prepare BAM files for merging per sample
        ch_bam_sample_to_merge = RM_TRANSITIONS.out.rm_trans_bam.map { meta, bam ->
            // update only the 'id' field in meta, keep all other fields
            [['id': meta.id.split("_")[0]], bam]
        }.groupTuple()
    } else { // if params.mapdamage2_rescale is false and params.remove_transitions is false, use deduplicated BAM files merged per library/PCR
        // Merge BAM files per sample
        ch_bam_sample_to_merge = SAMREMOVEDUP_LIB.out.dedup.map { meta, bam ->
            // update only the 'id' field in meta, keep all other fields
            [['id': meta.id.split("_")[0]], bam]
        }.groupTuple()
    }

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


    ////////////////////////////////////////////////////////////////////////////
    // GATK Indel Realignment
    ////////////////////////////////////////////////////////////////////////////

    if (params.indel_realignment.toBoolean()) {
        GATK_INDEL_REALIGNER(
            bam,
            reference,
        )
        ch_bam      = GATK_INDEL_REALIGNER.out.bam
        ch_versions = ch_versions.mix(GATK_INDEL_REALIGNER.out.versions)
    }

    emit:
    mq_filtered_bam                             = SAMTOOLS_VIEW_MQ.out.bam                                                                  // channel: [ val(meta), [ bam ] ]
    mq_filtered_index                           = SAMTOOLS_VIEW_MQ_INDEX.out.bai                                                            // channel: [ val(meta), [ bai ] ]
    rm_short_reads_bam                          = params.readlength == "auto" ? RM_SHORT_READS.out.bam : Channel.empty()                    // channel: [ val(meta), [ bam ] ]
    rm_short_reads_index                        = params.readlength == "auto" ? RM_SHORT_READS_INDEX.out.bai : Channel.empty()              // channel: [ val(meta), [ bai ] ]
    merged_bam_lib                              = SAMTOOLS_MERGE_LIB.out.bam                                                                // channel: [ val(meta), [ bam ] ]
    merged_bam_lib_index                        = SAMTOOLS_MERGE_LIB_INDEX.out.bai                                                          // channel: [ val(meta), [ bai ] ]
    dedup_lib                                   = SAMREMOVEDUP_LIB.out.dedup                                                                // channel: [ val(meta), [ bam ] ]
    dedup_lib_index                             = SAMREMOVEDUP_LIB_INDEX.out.bai                                                            // channel: [ val(meta), [ bai ] ]
    mapdamage2_rescaled_bam                     = params.mapdamage2_rescale.toBoolean() ? MAPDAMAGE2.out.rescaled_bam : Channel.empty()     // channel: [ val(meta), [ bam ] ]
    mapdamage2_rescaled_index                   = params.mapdamage2_rescale.toBoolean() ? MAPDAMAGE2_INDEX.out.bai : Channel.empty()        // channel: [ val(meta), [ bai ] ]
    mapdamage2_fragmisincorporation_plot        = MAPDAMAGE2.out.fragmisincorporation_plot                                                  // channel: [ val(meta), path(fragmisincorporation_plot) ]
    mapdamage2_length_plot                      = MAPDAMAGE2.out.length_plot                                                                // channel: [ val(meta), path(length_plot) ]
    mapdamage2_misincorporation                 = MAPDAMAGE2.out.misincorporation                                                           // channel: [ val(meta), path(misincorporation) ]
    mapdamage2_lgdistribution                   = MAPDAMAGE2.out.lgdistribution                                                             // channel: [ val(meta), path(lgdistribution) ]
    mapdamage2_dnacomp                          = MAPDAMAGE2.out.dnacomp                                                                    // channel: [ val(meta), path(dnacomp) ]
    mapdamage2_stats_out_mcmc_hist              = MAPDAMAGE2.out.stats_out_mcmc_hist                                                        // channel: [ val(meta), path(stats_out_mcmc_hist) ]
    mapdamage2_stats_out_mcmc_iter              = MAPDAMAGE2.out.stats_out_mcmc_iter                                                        // channel: [ val(meta), path(stats_out_mcmc_iter) ]
    mapdamage2_stats_out_mcmc_trace             = MAPDAMAGE2.out.stats_out_mcmc_trace                                                       // channel: [ val(meta), path(stats_out_mcmc_trace) ]
    mapdamage2_stats_out_mcmc_iter_summ_stat    = MAPDAMAGE2.out.stats_out_mcmc_iter_summ_stat                                              // channel: [ val(meta), path(stats_out_mcmc_iter_summ_stat) ]
    mapdamage2_stats_out_mcmc_post_pred         = MAPDAMAGE2.out.stats_out_mcmc_post_pred                                                   // channel: [ val(meta), path(stats_out_mcmc_post_pred) ]
    mapdamage2_stats_out_mcmc_correct_prob      = MAPDAMAGE2.out.stats_out_mcmc_correct_prob                                                // channel: [ val(meta), path(stats_out_mcmc_correct_prob) ]
    mapdamage2_dnacomp_genome                   = MAPDAMAGE2.out.dnacomp_genome                                                             // channel: [ val(meta), path(dnacomp_genome) ]
    mapdamage2_pctot_freq                       = MAPDAMAGE2.out.pctot_freq                                                                 // channel: [ val(meta), path(pctot_freq) ]
    mapdamage2_pgtoa_freq                       = MAPDAMAGE2.out.pgtoa_freq                                                                 // channel: [ val(meta), path(pgtoa_freq) ]
    mapdamage2_fasta                            = MAPDAMAGE2.out.fasta                                                                      // channel: [ val(meta), path(fasta) ]
    mapdamage2_folder                           = MAPDAMAGE2.out.folder                                                                     // channel: [ val(meta), path(folder) ]
    rm_trans_bam                                = params.remove_transitions.toBoolean() ? RM_TRANSITIONS.out.rm_trans_bam : Channel.empty() // channel: [ val(meta), [ bam ] ]
    rm_trans_index                              = params.remove_transitions.toBoolean() ? RM_TRANSITIONS_INDEX.out.bai : Channel.empty()    // channel: [ val(meta), [ bai ] ]
    merged_bam_sample                           = SAMTOOLS_MERGE_SAMPLE.out.bam                                                             // channel: [ val(meta), [ bam ] ]
    merged_bam_sample_index                     = SAMTOOLS_MERGE_SAMPLE_INDEX.out.bai                                                       // channel: [ val(meta), [ bai ] ]
    dedup_sample                                = SAMREMOVEDUP_SAMPLE.out.dedup                                                             // channel: [ val(meta), [ bam ] ]
    dedup_sample_index                          = SAMREMOVEDUP_SAMPLE_INDEX.out.bai                                                         // channel: [ val(meta), [ bai ] ]
    versions                                    = ch_versions                                                                               // channel: [ versions.yml ]
}