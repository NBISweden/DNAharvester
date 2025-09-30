#! /usr/bin/env nextflow

include { SAMTOOLS_UNMAPPED_READS                   } from '../../../modules/local/samtools/samtools_unmapped_reads.nf'
include { BWA_INDEX as MS_BWA_INDEX                 } from '../../../modules/local/bwa/index.nf'
include { BOWTIE2_BUILD as MS_BOWTIE2_BUILD         } from '../../../modules/local/bowtie2/bowtie2_build.nf'
include { BWA_ALN as MS_BWA_ALN                     } from '../../../modules/local/bwa/aln.nf'
include { BWA_SAMSE as MS_BWA_SAMSE                 } from '../../../modules/local/bwa/samse.nf'
include { BWA_ALN_MEM as MS_BWA_ALN_MEM             } from '../../../modules/local/bwa/bwa_aln_mem.nf'
include { BOWTIE2 as MS_BOWTIE2                     } from '../../../modules/local/bowtie2/bowtie2.nf'
include { SAMTOOLS_VIEW_MQ as MS_SAMTOOLS_VIEW_MQ   } from '../../../modules/local/samtools/samtools_view_mq.nf'
include { SAMTOOLS_MERGE as MS_SAMTOOLS_MERGE       } from '../../../modules/local/samtools/samtools_merge.nf'
include { SAMREMOVEDUP as MS_SAMREMOVEDUP           } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_INDEX as MS_SAMREMOVEDUP_INDEX   } from '../../../modules/nf-core/samtools/index/main'
include { FILTERBAM as MS_FILTERBAM                 } from '../../../modules/local/stats_output/filterbam.nf'
include { FILTERBAM_PLOT as MS_FILTERBAM_PLOT       } from '../../../modules/local/microbial_screening/filterbam_plot.nf'


workflow MICROBIAL_SCREENING {
    take:
    raw_bam
    ms_reference

    main:
    ch_versions             = Channel.empty()

    // get unmapped reads from the raw BAM file
    SAMTOOLS_UNMAPPED_READS ( raw_bam )
    ch_versions             = ch_versions.mix(SAMTOOLS_UNMAPPED_READS.out.versions)
    ch_unmapped_reads       = SAMTOOLS_UNMAPPED_READS.out.unmapped_fastq

    ////////////////////////////////////////////////////////////////////////////
    // Index the reference database if the index files are not present in the dir
    ////////////////////////////////////////////////////////////////////////////

    if (params.ms_mapping_tool == 'bwa-aln' || params.ms_mapping_tool == 'bwa-aln-mem') {
        // Build the BWA index
        MS_BWA_INDEX (ms_reference, file(params.ms_reference).getParent())
        ch_versions         = ch_versions.mix(MS_BWA_INDEX.out.versions)
        ch_ms_reference_index  = MS_BWA_INDEX.out.index_dir
    } else if (params.ms_mapping_tool == 'bowtie2') {
        // Build the Bowtie2 index
        MS_BOWTIE2_BUILD (ms_reference, file(params.ms_reference).getParent())
        ch_versions         = ch_versions.mix(MS_BOWTIE2_BUILD.out.versions)
        ch_ms_reference_index  = MS_BOWTIE2_BUILD.out.index_dir
    } else {
        error "Invalid mapping tool specified: ${params.ms_mapping_tool}. Use 'bwa-aln', 'bwa-aln-mem', or 'bowtie2'."
    }

    ////////////////////////////////////////////////////////////////////////////
    // Map the unmapped reads to the microbial reference database
    ////////////////////////////////////////////////////////////////////////////

    if (params.ms_mapping_tool == 'bwa-aln') {
        MS_BWA_ALN ( ch_unmapped_reads, ch_ms_reference_index )
        ch_versions         = ch_versions.mix(MS_BWA_ALN.out.versions)
        ch_bwa_samse        = ch_unmapped_reads.join(MS_BWA_ALN.out.sai)

        MS_BWA_SAMSE ( ch_bwa_samse, ch_ms_reference_index )
        ch_versions         = ch_versions.mix(MS_BWA_SAMSE.out.versions)
        ch_raw_bam          = MS_BWA_SAMSE.out.bam

    } else if (params.ms_mapping_tool == 'bwa-aln-mem') {
        MS_BWA_ALN_MEM ( ch_unmapped_reads, ch_ms_reference_index )
        ch_versions         = ch_versions.mix(MS_BWA_ALN_MEM.out.versions)
        ch_raw_bam          = MS_BWA_ALN_MEM.out.bam
    } else if (params.ms_mapping_tool == 'bowtie2') {
        MS_BOWTIE2 ( ch_unmapped_reads, ch_ms_reference_index )
        ch_versions         = ch_versions.mix(MS_BOWTIE2.out.versions)
        ch_raw_bam          = MS_BOWTIE2.out.bam
    } else {
        error "Invalid mapping tool specified: ${params.ms_mapping_tool}. Use 'bwa-aln' or 'bwa-aln-mem' or 'bowtie2'."
    }

    ////////////////////////////////////////////////////////////////////////////
    // Filter and Merge the BAM files
    ////////////////////////////////////////////////////////////////////////////

    // mapping quality filter
    MS_SAMTOOLS_VIEW_MQ ( ch_raw_bam )
    ch_versions             = ch_versions.mix ( MS_SAMTOOLS_VIEW_MQ.out.versions )

    // Prepare channel to merge BAM files per sample
    ch_ms_bams_per_sample  = MS_SAMTOOLS_VIEW_MQ.out.bam.map { meta, bam ->
        // update only the 'id' field in meta, keep all other fields
        [['id': meta.id.split("_")[0]], bam]
        }.groupTuple()

    // Merge BAM files per sample
    MS_SAMTOOLS_MERGE ( ch_ms_bams_per_sample, ms_reference )
    ch_versions             = ch_versions.mix(MS_SAMTOOLS_MERGE.out.versions)

    // Remove PCR duplicates
    MS_SAMREMOVEDUP ( MS_SAMTOOLS_MERGE.out.bam, ms_reference )
    ch_versions             = ch_versions.mix(MS_SAMREMOVEDUP.out.versions)
    MS_SAMREMOVEDUP_INDEX ( MS_SAMREMOVEDUP.out.dedup )
    ch_versions             = ch_versions.mix(MS_SAMREMOVEDUP_INDEX.out.versions)
    ch_bam_bai              = MS_SAMREMOVEDUP.out.dedup.join(MS_SAMREMOVEDUP_INDEX.out.bai)

    ////////////////////////////////////////////////////////////////////////////
    // Generate stats and plots
    ////////////////////////////////////////////////////////////////////////////

    // run filterBAM to generate stats
    MS_FILTERBAM ( ch_bam_bai )
    ch_versions             = ch_versions.mix ( MS_FILTERBAM.out.versions )

    // generate plot from filterBAM output
    MS_FILTERBAM_PLOT ( ch_bam_bai, MS_FILTERBAM.out.filterBAM_stats )
    ch_versions             = ch_versions.mix ( MS_FILTERBAM_PLOT.out.versions )


    emit:
    unmapped_fastq                = SAMTOOLS_UNMAPPED_READS.out.unmapped_fastq  // channel: [ val(meta), [ fastq ] ]
    ms_reference_index            = ch_ms_reference_index                       // channel: path(index)
    ms_raw_bam                    = ch_raw_bam                                  // channel: [ val(meta), [ bam ] ]
    ms_bam_bai                    = ch_bam_bai                                  // channel: [ val(meta), [ bai ] ]
    ms_filterBAM_stats            = MS_FILTERBAM.out.filterBAM_stats            // channel: [ val(meta), [ filterBAM.csv ] ]
    ms_screening_plot             = MS_FILTERBAM_PLOT.out.ms_plot               // channel: [ val(meta), [ pathogen_screening_plot.pdf ] ]
    versions                      = ch_versions                                 // channel: [ versions.yml ]
}