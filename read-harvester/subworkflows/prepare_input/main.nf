#! /usr/bin/env nextflow

// Check input samplesheet and get read channels
// Modified from nf-core pipeline template 
// https://github.com/nf-core/tools/blob/e5ce6ce20304835bd40f102f038b7e1aadc888b2/nf_core/pipeline-template/subworkflows/local/input_check.nf

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'
//include { READ_YAML } from '../../modules/local/read_yaml'

workflow PREPARE_INPUT {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    //paramsyaml

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

//    READ_YAML ( paramsyaml )
// add code to read the params.yaml file with e.g. info about reference genome path

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()

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