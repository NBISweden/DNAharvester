#! /usr/bin/env nextflow

// Analysis script for Read harvester pipeline

nextflow.enable.dsl = 2

// Import subworkflows
include { INPUT_CHECK                } from "$projectDir/subworkflows/local/input_check/main"
include { FASTQ_PROCESSING           } from "$projectDir/subworkflows/local/fastq_processing/main"
include { PROCESSED_FASTQ_QC         } from "$projectDir/subworkflows/local/processed_fastq_qc/main"
include { MAPPING                    } from "$projectDir/subworkflows/local/mapping/main"
include { COMPETITIVE_MAPPING        } from "$projectDir/subworkflows/local/competitive_mapping/main"
include { RAW_BAM_QC                 } from "$projectDir/subworkflows/local/raw_bam_qc/main"
include { BAM_PROCESSING             } from "$projectDir/subworkflows/local/bam_processing/main"
include { PROCESSED_BAM_QC           } from "$projectDir/subworkflows/local/processed_bam_qc/main"
include { RANDOM_SAMPLING_BAM        } from "$projectDir/subworkflows/local/random_sampling_bam/main"
include { STATS_OUTPUT               } from "$projectDir/subworkflows/local/stats_output/main"

workflow {

    // Define workflow stages
    def recognized_workflow_stages = ['fastq_processing','mapping','processed_fastq_qc','raw_bam_qc', 'bam_processing', 'processed_bam_qc', 'random_sampling_bam','output_stats']

    // Check input
    def workflow_steps = params.steps.tokenize(",")
    if ( ! workflow_steps.every { it in recognized_workflow_stages } ) {
        error "Unrecognised workflow step in $params.steps ( $recognized_workflow_stages )"
    }

    // The primary workflow for the read harvester pipeline
    log.info("""
    Running Read Harvester.
    """)

    // Read in data and create channels
    INPUT_CHECK ( params.samplesheet )

    ch_reference = Channel.fromPath( params.reference, checkIfExists: true )
        .map { it -> [[id:it.Name], it] }.collect()
    
    ch_competitive_reference = params.competitive_reference ? Channel.fromPath( params.competitive_reference, checkIfExists: true )
        .map { it -> [[id:it.Name], it] }.collect() : Channel.empty()

    ch_competitive_reference_index = params.competitive_reference ? Channel.fromFilePairs("${params.competitive_reference}*.{amb,ann,bwt,pac,sa}", size: 5, checkIfExists: true)
        .map { id, files ->
            def parentDir = files[0].getParent()
            return [[id:id], parentDir] } : Channel.empty()

    ch_intervals = params.intervals ? Channel.fromPath( params.intervals, checkIfExists: true )
        .map { it -> [[id:it.Name], it] }.collect() : Channel.empty()

    // Merge paired-end reads, trim adapters and filter for minimum read length
    if ( 'fastq_processing' in workflow_steps ) {
        FASTQ_PROCESSING (
            INPUT_CHECK.out.reads
        )
    }

    // Run FastQC, MultiQC and read statistics on processed reads
    if ( 'processed_fastq_qc' in workflow_steps ) {
        PROCESSED_FASTQ_QC (
            FASTQ_PROCESSING.out.reads,
            FASTQ_PROCESSING.out.json
        )
    }

    // Map with bwa-aln (aDNA parameters) and convert to bam
    if ( 'mapping' in workflow_steps ) {
        // Competitive mapping to a concatenated reference (target plus decoy)
        if (params.competitive_reference && file( params.competitive_reference ).exists()) {
            COMPETITIVE_MAPPING (
                    ch_competitive_reference,
                    ch_competitive_reference_index,
                    ch_reference,
                    FASTQ_PROCESSING.out.reads
            )
        // Map to the reference genome assembly
        } else {
            MAPPING (
                ch_reference,
                FASTQ_PROCESSING.out.reads
            )
        }
    }

    // Run samtools flagstat, MapDamage2, AMBER and MultiQC on raw bam files
    if ( 'raw_bam_qc' in workflow_steps ) {
        RAW_BAM_QC (
            ch_reference,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.bam : MAPPING.out.bam,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.bai : MAPPING.out.bai,
        )
    }
    // Merge bam files per index, remove duplicates, merge bam files per sample, remove duplicates, realign indels
    if ( 'bam_processing' in workflow_steps ) {
        BAM_PROCESSING (
            ch_reference,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.fai : MAPPING.out.fai,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.bam : MAPPING.out.bam,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.bai : MAPPING.out.bai,
            RAW_BAM_QC.out.amber_txt
        )
    }

    // Run flagstat and MultiQC on processed bam files
    if ( 'processed_bam_qc' in workflow_steps ) {
        PROCESSED_BAM_QC (
            ch_reference,
            BAM_PROCESSING.out.mq_filtered_bam,
            BAM_PROCESSING.out.mq_filtered_index,
            BAM_PROCESSING.out.rm_short_reads_bam,
            BAM_PROCESSING.out.rm_short_reads_index,
            BAM_PROCESSING.out.merged_bam_lib,
            BAM_PROCESSING.out.merged_bam_lib_index,
            BAM_PROCESSING.out.dedup_lib,
            BAM_PROCESSING.out.dedup_lib_index,
            BAM_PROCESSING.out.merged_bam_sample,
            BAM_PROCESSING.out.merged_bam_sample_index,
            BAM_PROCESSING.out.dedup_sample,
            BAM_PROCESSING.out.dedup_sample_index,
            BAM_PROCESSING.out.realigned,
            ch_intervals
        )
    }

    // Run ANGSD -doHaploCall 1 to sample a random base at each site from bam files
    if ( 'random_sampling_bam' in workflow_steps ) {
        RANDOM_SAMPLING_BAM (
            BAM_PROCESSING.out.realigned,
            ch_reference,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.fai : MAPPING.out.fai
        )
    }

    // Output stats
    if ( 'output_stats' in workflow_steps ) {
        STATS_OUTPUT (
            INPUT_CHECK.out.reads,
            FASTQ_PROCESSING.out.fastp_log,
            RAW_BAM_QC.out.flagstat,
            PROCESSED_BAM_QC.out.dedup_lib_flagstat,
            BAM_PROCESSING.out.realigned
        )
    }

}

workflow.onComplete {
    if( workflow.success ){
        log.info("""
        Thank you for using Read Harvester.

        Results are located in the folder: $params.outdir
        """)
    } else {
        log.info("""
        The pipeline completed unsuccessfully.

        Please read the error message. If you need help to solve your issue,
        feel free to reach out via slack or by opening an issue at
        https://github.com/NBISweden/LTS-L_Dalen_2302/issues.
        """)
    }
}