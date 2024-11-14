#! /usr/bin/env nextflow

// Analysis script for Read harvester pipeline

nextflow.enable.dsl = 2

// Import subworkflows
include { INPUT_CHECK                } from "$projectDir/subworkflows/local/input_check/main"
include { FASTQ_PROCESSING           } from "$projectDir/subworkflows/local/fastq_processing/main"
include { MAPPING                    } from "$projectDir/subworkflows/local/mapping/main"
include { PROCESSED_FASTQ_RAW_BAM_QC } from "$projectDir/subworkflows/local/processed_fastq_raw_bam_qc/main"
include { BAM_PROCESSING             } from "$projectDir/subworkflows/local/bam_processing/main"
include { PROCESSED_BAM_QC           } from "$projectDir/subworkflows/local/processed_bam_qc/main"
include { RANDOM_SAMPLING_BAM        } from "$projectDir/subworkflows/local/random_sampling_bam/main"

workflow {

    // Define workflow stages
    def recognized_workflow_stages = ['fastq_processing','mapping','processed_fastq_raw_bam_qc', 'bam_processing', 'processed_bam_qc', 'random_sampling_bam']

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
    Channel.fromPath( params.reference, checkIfExists: true )
        .set{ reference }

    // Merge paired-end reads, trim adapters and filter for minimum read length
    if ( 'fastq_processing' in workflow_steps ) {
        FASTQ_PROCESSING (
            INPUT_CHECK.out.reads
        )
    }

    // Index the reference genome, map with bwa-aln (aDNA parameters) and convert to bam
    if ( 'mapping' in workflow_steps ) {
        MAPPING (
            params.reference ? file( params.reference, checkIfExists: true ) : [],
            FASTQ_PROCESSING.out.reads
        ) 
    }

    // Run FastQC, samtools flagstat, MapDamage2, AMBER and MultiQC to assess the data quality
    if ( 'processed_fastq_raw_bam_qc' in workflow_steps ) {
        PROCESSED_FASTQ_RAW_BAM_QC (
            params.reference ? file( params.reference, checkIfExists: true ) : [],
            FASTQ_PROCESSING.out.reads,
            FASTQ_PROCESSING.out.json,
            MAPPING.out.bam,
            MAPPING.out.bai
        )
    }

    // Index the reference genome, merge bam files per index, remove duplicates, merge bam files per sample, remove duplicates, realign indels
    if ( 'bam_processing' in workflow_steps ) {
        BAM_PROCESSING (
            params.reference ? file( params.reference, checkIfExists: true ) : [],
            MAPPING.out.bam
        ) 
    }

    // Run QualiMap and MultiQC on processed bam files
    if ( 'processed_bam_qc' in workflow_steps ) {
        PROCESSED_BAM_QC (
            params.reference ? file( params.reference, checkIfExists: true ) : [],
            BAM_PROCESSING.out.merged_bam_lib,
            BAM_PROCESSING.out.merged_bam_lib_index,
            BAM_PROCESSING.out.dedup_lib,
            BAM_PROCESSING.out.dedup_lib_index,
            BAM_PROCESSING.out.merged_bam_sample,
            BAM_PROCESSING.out.merged_bam_sample_index,
            BAM_PROCESSING.out.dedup_sample,
            BAM_PROCESSING.out.dedup_sample_index,
            BAM_PROCESSING.out.realigned,
            params.intervals
        ) 
    }

    // Run ANGSD -doHaploCall 1 to sample a random base at each site from bam files
    if ( 'random_sampling_bam' in workflow_steps ) {
        RANDOM_SAMPLING_BAM (
            BAM_PROCESSING.out.realigned,
            params.reference ? file( params.reference, checkIfExists: true ) : [],
            BAM_PROCESSING.out.fai
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