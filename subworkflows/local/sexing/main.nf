#! /usr/bin/env nextflow

/*
 * SEXING SUBWORKFLOW
 *
 * Determines biological sex from processed BAM files using:
 * 1. X-to-autosome read ratio method (default)
 * 2. Y-chromosome read presence method (optional)
 *
 * Required parameters:
 *   - params.sexing_x_chr: X chromosome name(s), space-separated
 *   - params.sexing_autosomes: Autosome name(s) for normalization, space-separated
 *
 * Optional parameters:
 *   - params.sexing_y_chr: Y chromosome name (enables Y-chr method)
 *   - params.sexing_method: 'x_ratio', 'y_chr', or 'both' (default: 'x_ratio')
 *   - params.sexing_x_ratio_threshold: Female classification threshold (default: 0.75)
 *   - params.sexing_min_reads: Minimum reads for reliable determination (default: 1000)
 */

include { SAMTOOLS_IDXSTATS as S_SAMTOOLS_IDXSTATS } from '../../../modules/local/samtools/samtools_idxstats.nf'
include { SEX_DETERMINATION                        } from '../../../modules/local/sexing/sex_determination.nf'
include { MERGE_SEX_REPORTS                        } from '../../../modules/local/sexing/merge_sex_reports.nf'


workflow SEXING {
    take:
    processed_bam

    main:
    ch_versions = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Run samtools idxstats to get the number of reads mapped to each reference sequence
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    S_SAMTOOLS_IDXSTATS(processed_bam)
    ch_versions = ch_versions.mix(S_SAMTOOLS_IDXSTATS.out.versions)
    ch_idxstats = S_SAMTOOLS_IDXSTATS.out.idxstats

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. Calculate sex determination from idxstats
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    SEX_DETERMINATION(
        ch_idxstats,
        params.sexing_x_chr ?: '',
        params.sexing_autosomes ?: '',
        params.sexing_y_chr ?: ''
    )
    ch_versions = ch_versions.mix(SEX_DETERMINATION.out.versions)
    ch_sex_reports = SEX_DETERMINATION.out.sex_report

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 3. Merge all individual sex reports into a summary
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    MERGE_SEX_REPORTS(
        ch_sex_reports.map { meta, report -> report }.collect()
    )
    ch_versions = ch_versions.mix(MERGE_SEX_REPORTS.out.versions)

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    sex_reports     = ch_sex_reports                    // channel: [ meta, sex_determination.tsv ]
    summary         = MERGE_SEX_REPORTS.out.summary     // channel: [ sex_determination_summary.tsv ]
    versions        = ch_versions                       // channel: [ versions.yml ]
}