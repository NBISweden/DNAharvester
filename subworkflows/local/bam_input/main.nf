#! /usr/bin/env nextflow

include { SAMTOOLS_FAIDX    as BI_SAMTOOLS_FAIDX    } from '../../../modules/local/samtools/samtools_faidx.nf'
include { SAMTOOLS_CHECK_SQ as BI_SAMTOOLS_CHECK_SQ } from '../../../modules/local/samtools/samtools_check_sq.nf'
include { SAMTOOLS_INDEX    as BI_SAMTOOLS_INDEX    } from '../../../modules/nf-core/samtools/index/main'


workflow BAM_INPUT {
    take:
    bam         // channel: [ val(meta), path(bam), path(bai) or [] ]
    reference

    main:
    ch_versions = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Index the reference genome and check that the BAM files were mapped to it
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    BI_SAMTOOLS_FAIDX ( reference )
    ch_versions = ch_versions.mix(BI_SAMTOOLS_FAIDX.out.versions)

    BI_SAMTOOLS_CHECK_SQ (
        bam.map { meta, bam_f, bai -> [meta, bam_f] },
        BI_SAMTOOLS_FAIDX.out.fai
    )
    ch_versions = ch_versions.mix(BI_SAMTOOLS_CHECK_SQ.out.versions)

    // ANGSD needs identical @SQ headers in all BAM files for joint variant calling.
    // The BAM files are passed on only after this check, so no tasks are running when it exits
    ch_sq_checked = BI_SAMTOOLS_CHECK_SQ.out.sq
        .map { meta, sq -> [meta.id, sq.text] }
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { sq_list ->
            def differing = sq_list.findAll { it[1] != sq_list[0][1] }.collect { it[0] }
            if ( differing && params.variant_calling.toBoolean() && params.variant_calling_angsd.toBoolean() ) {
                log.error """`variant_calling_angsd` requires identical @SQ headers in all BAM files, but the @SQ header of
                ${differing.join(', ')} differs from ${sq_list[0][0]} (e.g. BAM files mapped with and without a competitive reference genome).
                Please provide BAM files with identical @SQ headers or use `variant_calling_bcftools` instead.
                Exiting the pipeline......!
                """
                System.exit(1)
            }
            true
        }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. Index the BAM files without a BAI file next to them
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ch_bam_branched = BI_SAMTOOLS_CHECK_SQ.out.bam
        .combine(ch_sq_checked)
        .map { meta, bam_f, sq_checked -> [meta, bam_f] }
        .join(bam.map { meta, bam_f, bai -> [meta, bai] })
        .branch { meta, bam_f, bai ->
            indexed:     bai
            not_indexed: true
        }

    BI_SAMTOOLS_INDEX ( ch_bam_branched.not_indexed.map { meta, bam_f, bai -> [meta, bam_f] } )
    ch_versions = ch_versions.mix(BI_SAMTOOLS_INDEX.out.versions)

    ch_bam_bai = ch_bam_branched.indexed
        .mix(ch_bam_branched.not_indexed.map { meta, bam_f, bai -> [meta, bam_f] }.join(BI_SAMTOOLS_INDEX.out.bai))

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    bam         = ch_bam_bai.map { meta, bam_f, bai -> [meta, bam_f] }  // channel: [ val(meta), bam ]
    bai         = ch_bam_bai.map { meta, bam_f, bai -> [meta, bai] }    // channel: [ val(meta), bai ]
    fai         = BI_SAMTOOLS_FAIDX.out.fai                             // channel: [ val(meta), fai ]
    versions    = ch_versions                                           // channel: [ versions.yml ]
}
