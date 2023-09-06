#! /usr/bin/env nextflow

// Analysis script for Read harvester pipeline

nextflow.enable.dsl = 2

include { INPUT_CHECK } from "$projectDir/subworkflows/input_check/main"
//include { MERGE_FILTER_READS   } from "$projectDir/subworkflows/merge_filter_reads/main"


workflow {

    // Define workflow stages
    def recognized_workflow_stages = ['read_processing']

    // Check input
    def workflow_steps = params.steps.tokenize(",")
    if ( ! workflow_steps.every { it in recognized_workflow_stages } ) {
        error "Unrecognised workflow step in $params.steps ( $recognized_workflow_stages )"
    }

    // The primary workflow for the read harvester pipeline
    log.info("""
    Running Read Harvester.
    """)

    // Read in data
    INPUT_CHECK ( params.samplesheet )

    // Merge paired-end reads, trim adapters and filter for minimum read length
    //if ( 'read_processing' in workflow_steps ) {
    //    FASTP (
    //        INPUT_CHECK.out // specify channel
    //    ) 
    //}

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