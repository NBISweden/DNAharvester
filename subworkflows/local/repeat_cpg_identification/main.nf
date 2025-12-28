#! /usr/bin/env nextflow

include { REPEATMODELER         as RCI_REPEATMODELER         } from '../../../modules/local/repeat_cpg_identification/repeatmodeler'
include { REPEATMASKER          as RCI_REPEATMASKER          } from '../../../modules/local/repeat_cpg_identification/repeatmasker'
include { SAMTOOLS_FAIDX        as RCI_SAMTOOLS_FAIDX        } from '../../../modules/local/samtools/samtools_faidx.nf'
include { CREATE_REPEATS_BED    as RCI_CREATE_REPEATS_BED    } from '../../../modules/local/repeat_cpg_identification/create_repeats_bed'
include { CREATE_REPMA_BED      as RCI_CREATE_REPMA_BED      } from '../../../modules/local/repeat_cpg_identification/create_repma_bed'
include { CREATE_CPG_BED        as RCI_CREATE_CPG_BED        } from '../../../modules/local/repeat_cpg_identification/create_cpg_bed'
include { CREATE_REPMA_CPG_BED  as RCI_CREATE_REPMA_CPG_BED  } from '../../../modules/local/repeat_cpg_identification/create_repma_cpg_bed'


workflow REPEAT_CPG_IDENTIFICATION {
    take:
    reference

    main:
    ch_versions = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Identify repeats in the reference genome
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Run RepeatModeler
    RCI_REPEATMODELER ( reference )
    ch_versions = ch_versions.mix(RCI_REPEATMODELER.out.versions)

    // Run RepeatMasker
    RCI_REPEATMASKER ( reference, RCI_REPEATMODELER.out.upper_fasta, RCI_REPEATMODELER.out.consensi)
    ch_versions = ch_versions.mix(RCI_REPEATMASKER.out.versions)

    // Create the repeat bed file
    RCI_CREATE_REPEATS_BED ( reference, RCI_REPEATMASKER.out.repeatmasker_out )
    ch_versions = ch_versions.mix(RCI_CREATE_REPEATS_BED.out.versions)

    // Index the reference genome
    RCI_SAMTOOLS_FAIDX ( reference )
    ch_versions = ch_versions.mix(RCI_SAMTOOLS_FAIDX.out.versions)

    // Create the repeats masked (repma) bed file
    RCI_CREATE_REPMA_BED ( reference, RCI_SAMTOOLS_FAIDX.out.fai, RCI_CREATE_REPEATS_BED.out.repeats_bed )
    ch_versions = ch_versions.mix(RCI_CREATE_REPMA_BED.out.versions)

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. Identify CpG sites and create combined CpG and repeat masked files
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create the CpG bed file
    RCI_CREATE_CPG_BED ( reference )
    ch_versions = ch_versions.mix(RCI_CREATE_CPG_BED.out.versions)

    // Merge the CpG and repeat bed files
    RCI_CREATE_REPMA_CPG_BED (
        reference,
        RCI_CREATE_REPMA_BED.out.genomefile,
        RCI_CREATE_REPMA_BED.out.sorted_repeats_bed,
        RCI_CREATE_REPMA_BED.out.ref_bed,
        RCI_CREATE_CPG_BED.out.cpg_bed)
    ch_versions = ch_versions.mix(RCI_CREATE_REPMA_CPG_BED.out.versions)

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    sorted_repeats_bed      = RCI_CREATE_REPMA_BED.out.sorted_repeats_bed           // channel: path(sorted_repeats_bed)
    repma_bed               = RCI_CREATE_REPMA_BED.out.repma_bed                    // channel: path(repma_bed)
    cpg_bed                 = RCI_CREATE_CPG_BED.out.cpg_bed                        // channel: path(cpg_bed)
    no_cpg_bed              = RCI_CREATE_REPMA_CPG_BED.out.no_cpg_bed               // channel: path(noCpG_bed)
    cpg_repeats_sorted_bed  = RCI_CREATE_REPMA_CPG_BED.out.cpg_repeats_sorted_bed   // channel: path(cpg_repeats_sorted_bed)
    no_cpg_repma_bed        = RCI_CREATE_REPMA_CPG_BED.out.no_cpg_repma_bed         // channel: path(noCpG_repma_bed)
    versions                = ch_versions                                           // channel: [ versions.yml ]
}