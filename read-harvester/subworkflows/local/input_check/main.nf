#! /usr/bin/env nextflow

// Check input samplesheet and get read channels
// Modified from nf-core pipeline template 
// https://github.com/nf-core/tools/blob/e5ce6ce20304835bd40f102f038b7e1aadc888b2/nf_core/pipeline-template/subworkflows/local/input_check.nf

include { SAMPLESHEET_CHECK  } from '../../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }
    

    emit:
    reads                                         // channel: [ val(meta), [ reads ] ]
    csv      = SAMPLESHEET_CHECK.out.csv          // channel: [ samplesheet.valid.csv ]
    versions = SAMPLESHEET_CHECK.out.versions     // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map: [id, single_end]
    // TO DO Update to take strandedness and damage treatment
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()
    // readgroup: ID = ID = readgroup id (flowcell-id.lane-nr.library-index-nr), SM = sample-id, PL = sequencing platform (e.g. Illumina, NovaSeq), LB = library-index-nr based on the number of unique libraries per sample in samplesheet
    meta.read_group    = "@RG\\tID:" + row.flowcell_id + "." + row.lane + "." + row.sample.split('_')[1] + "\\tSM:" + row.sample.split('_')[0] + "\\tPL:" + row.seq_platform + "\\tLB:" + row.sample.split('_')[1]

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}