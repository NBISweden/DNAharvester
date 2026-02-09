#! /usr/bin/env nextflow

/*
 * SEXING SUBWORKFLOW
 *
 * Determines biological sex from processed BAM files using:
 * 1. X-to-autosome read ratio method (when sexing_x_chr is provided)
 * 2. Y-chromosome read presence method (when sexing_y_chr is provided)
 *
 * Required parameters:
 *   - params.sexing_x_chr: X chromosome name(s), space-separated (for X-ratio method)
 *   - params.sexing_autosomes: Autosome name(s) for normalization, space-separated
 *
 * Optional parameters:
 *   - params.sexing_y_chr: Y chromosome name(s) (enables Y-chr method)
 */

include { SAMTOOLS_IDXSTATS as S_SAMTOOLS_IDXSTATS  } from '../../../modules/local/samtools/samtools_idxstats.nf'
include { X_CHR_SEXING      as S_X_CHR_SEXING       } from '../../../modules/local/sexing/x_chr_sexing.nf'
include { Y_CHR_SEXING      as S_Y_CHR_SEXING       } from '../../../modules/local/sexing/y_chr_sexing.nf'
include { MERGE_SEX_REPORTS as S_MERGE_X_CHR_REPORTS  } from '../../../modules/local/sexing/merge_sex_reports.nf'
include { MERGE_SEX_REPORTS as S_MERGE_Y_CHR_REPORTS  } from '../../../modules/local/sexing/merge_sex_reports.nf'


workflow SEXING {
    take:
    processed_bam

    main:
    ch_versions = Channel.empty()
    ch_x_sex_reports = Channel.empty()
    ch_x_sex_summary = Channel.empty()
    ch_y_sex_reports = Channel.empty()
    ch_y_sex_summary = Channel.empty()
    ch_x_merged_summary = Channel.empty()
    ch_y_merged_summary = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Run samtools idxstats to get the number of reads mapped to each reference sequence
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    S_SAMTOOLS_IDXSTATS ( processed_bam )
    ch_versions = ch_versions.mix(S_SAMTOOLS_IDXSTATS.out.versions)
    ch_idxstats = S_SAMTOOLS_IDXSTATS.out.idxstats

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. X-chromosome sex determination (if sexing_x_chr is provided)
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if ( params.sexing_x_chr && params.sexing_autosomes ) {
        S_X_CHR_SEXING (
            ch_idxstats,
            params.sexing_x_chr,
            params.sexing_autosomes
        )
        ch_versions = ch_versions.mix(S_X_CHR_SEXING.out.versions)
        ch_x_sex_reports = S_X_CHR_SEXING.out.sex_report
        ch_x_sex_summary = S_X_CHR_SEXING.out.sex_summary

        // Merge X-chr sex reports
        S_MERGE_X_CHR_REPORTS (
            ch_x_sex_summary.map { meta, report -> report }.collect()
        )
        ch_versions = ch_versions.mix(MERGE_X_CHR_REPORTS.out.versions)
        ch_x_merged_summary = MERGE_X_CHR_REPORTS.out.summary
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 3. Y-chromosome sex determination (if sexing_y_chr is provided)
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // if ( params.sexing_y_chr && params.sexing_autosomes ) {
    //     Y_CHR_SEX_DETERMINATION (
    //         ch_idxstats,
    //         params.sexing_y_chr,
    //         params.sexing_autosomes
    //     )
    //     ch_versions = ch_versions.mix(Y_CHR_SEX_DETERMINATION.out.versions)
    //     ch_y_sex_reports = Y_CHR_SEX_DETERMINATION.out.sex_report
    //     ch_y_sex_summary = Y_CHR_SEX_DETERMINATION.out.sex_summary

    //     // Merge Y-chr sex reports
    //     MERGE_Y_CHR_REPORTS (
    //         ch_y_sex_summary.map { meta, report -> report }.collect()
    //     )
    //     ch_versions = ch_versions.mix(S_MERGE_Y_CHR_REPORTS.out.versions)
    //     ch_y_merged_summary = S_MERGE_Y_CHR_REPORTS.out.summary
    // }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    x_sex_reports       = ch_x_sex_reports          // channel: [ meta, x_chr_sexing.tsv ]
    x_sex_summary       = ch_x_sex_summary          // channel: [ meta, x_chr_sexing_summary.tsv ]
    x_merged_summary    = ch_x_merged_summary       // channel: [ sex_determination_summary.tsv ]
    y_sex_reports       = ch_y_sex_reports          // channel: [ meta, y_chr_sexing.tsv ]
    y_sex_summary       = ch_y_sex_summary          // channel: [ meta, y_chr_sexing_summary.tsv ]
    y_merged_summary    = ch_y_merged_summary       // channel: [ sex_determination_summary.tsv ]
    versions            = ch_versions               // channel: [ versions.yml ]
}