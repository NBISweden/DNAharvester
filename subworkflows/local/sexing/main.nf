#! /usr/bin/env nextflow

include { SAMTOOLS_IDXSTATS as S_SAMTOOLS_IDXSTATS  } from '../../../modules/local/samtools/samtools_idxstats.nf'
include { X_CHR_SEXING      as S_X_CHR_SEXING       } from '../../../modules/local/sexing/x_chr_sexing.nf'
include { Y_CHR_SEXING      as S_Y_CHR_SEXING       } from '../../../modules/local/sexing/y_chr_sexing.nf'


workflow SEXING {
    take:
    workflow_name
    processed_bam

    main:
    ch_versions = Channel.empty()

    // Initialize output channels as empty
    ch_x_chr_sexing_report      = Channel.empty()
    ch_x_chr_sexing_ploidy_plot = Channel.empty()
    ch_x_chr_sexing_summary     = Channel.empty()
    ch_y_chr_sexing_report      = Channel.empty()
    ch_merged_x_chr_sexing      = Channel.empty()
    ch_merged_y_chr_sexing      = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Run samtools idxstats to get the number of reads mapped to each reference sequence
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    S_SAMTOOLS_IDXSTATS ( processed_bam )
    ch_versions = ch_versions.mix(S_SAMTOOLS_IDXSTATS.out.versions)
    ch_idxstats = S_SAMTOOLS_IDXSTATS.out.idxstats

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. X-chromosome sex determination (if sexing_x_chr and sexing_autosomes are provided)
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if ( params.sexing_x_chr && params.sexing_autosomes ) {
        S_X_CHR_SEXING (
            ch_idxstats,
            params.sexing_x_chr,
            params.sexing_autosomes
        )
        ch_versions = ch_versions.mix(S_X_CHR_SEXING.out.versions)

        ch_x_chr_sexing_report      = S_X_CHR_SEXING.out.sexing_report
        ch_x_chr_sexing_ploidy_plot = S_X_CHR_SEXING.out.sexing_ploidy_plot
        ch_x_chr_sexing_summary     = S_X_CHR_SEXING.out.sexing_summary

        // Concatenate all output files
        ch_merged_x_chr_sexing = S_X_CHR_SEXING.out.sexing_summary
            .map { it[1] }
            .collectFile(name: "${workflow_name}_x_chr_sexing_summary.tsv", keepHeader: true, skip: 1, sort: true)
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 3. Y-chromosome sex determination (if sexing_y_chr and sexing_x_chr are provided)
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if ( params.sexing_y_chr && params.sexing_x_chr ) {
        S_Y_CHR_SEXING (
            ch_idxstats,
            params.sexing_y_chr,
            params.sexing_x_chr
        )
        ch_versions = ch_versions.mix(S_Y_CHR_SEXING.out.versions)

        ch_y_chr_sexing_report = S_Y_CHR_SEXING.out.sexing_report

        // Concatenate all output files
        ch_merged_y_chr_sexing = S_Y_CHR_SEXING.out.sexing_report
            .map { it[1] }
            .collectFile(name: "${workflow_name}_y_chr_sexing_summary.tsv", keepHeader: true, skip: 1, sort: true)
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    x_chr_sexing_report             = ch_x_chr_sexing_report
    x_chr_sexing_ploidy_plot        = ch_x_chr_sexing_ploidy_plot
    x_chr_sexing_summary            = ch_x_chr_sexing_summary
    y_chr_sexing_report             = ch_y_chr_sexing_report
    merged_x_chr_sexing             = ch_merged_x_chr_sexing
    merged_y_chr_sexing             = ch_merged_y_chr_sexing
    versions                        = ch_versions
}