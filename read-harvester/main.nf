#! /usr/bin/env nextflow

// Analysis script for Read harvester pipeline

nextflow.enable.dsl = 2

// Import subworkflows
include { INPUT_CHECK        } from "$projectDir/subworkflows/local/input_check/main"
include { MERGE_FILTER_READS } from "$projectDir/subworkflows/local/merge_filter_reads/main"
include { MAPPING            } from "$projectDir/subworkflows/local/mapping/main"
include { DATA_QC            } from "$projectDir/subworkflows/local/data_qc/main"

workflow {

    // Define workflow stages
    def recognized_workflow_stages = ['read_processing','mapping','data_qc']

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
    if ( 'read_processing' in workflow_steps ) {
        MERGE_FILTER_READS (
            INPUT_CHECK.out.reads
        )
    }

    // Index the reference genome, map with bwa-aln (aDNA parameters) and convert to bam
    if ( 'mapping' in workflow_steps ) {
        MAPPING (
            params.reference ? file( params.reference, checkIfExists: true ) : [],
            MERGE_FILTER_READS.out.reads
        ) 
    }

    // Run FastQC, MapDamage2, AMBER and MultiQC to assess the data quality
    if ( 'data_qc' in workflow_steps ) {
        DATA_QC (
            params.reference ? file( params.reference, checkIfExists: true ) : [],
            INPUT_CHECK.out.reads,
            MERGE_FILTER_READS.out.reads,
            MERGE_FILTER_READS.out.json,
            MAPPING.out.bam
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