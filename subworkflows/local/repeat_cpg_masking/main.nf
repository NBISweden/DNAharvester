#! /usr/bin/env nextflow

include { REPEATMODELER         } from '../../../modules/local/repeat_cpg_masking/repeatmodeler'
include { REPEATMASKER          } from '../../../modules/local/repeat_cpg_masking/repeatmasker'
include { SAMTOOLS_FAIDX        } from '../../../modules/local/samtools/faidx/main'
include { CREATE_REPEATS_BED    } from '../../../modules/local/repeat_cpg_masking/create_repeats_bed'
include { CREATE_REPMA_BED      } from '../../../modules/local/repeat_cpg_masking/create_repma_bed'
include { CREATE_CPG_BED        } from '../../../modules/local/repeat_cpg_masking/create_cpg_bed'
include { CREATE_REPMA_CPG_BED  } from '../../../modules/local/repeat_cpg_masking/create_repma_cpg_bed'


workflow REPEAT_CPG_MASKING {
    take:
    reference

    main:
    ch_versions = Channel.empty()

    // Run RepeatModeler
    REPEATMODELER ( reference )
    ch_versions = ch_versions.mix(REPEATMODELER.out.versions)

    // Run RepeatMasker
    REPEATMASKER ( REPEATMODELER.out.upper_fasta, REPEATMODELER.out.consensi)
    ch_versions = ch_versions.mix(REPEATMASKER.out.versions)

    // Index the reference genome
    SAMTOOLS_FAIDX ( reference )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // Create the repeat bed file
    CREATE_REPEATS_BED ( REPEATMASKER.out.repeatmasker_out )
    ch_versions = ch_versions.mix(CREATE_REPEATS_BED.out.versions)

    // Create the repeats masked (repma) bed file
    CREATE_REPMA_BED ( SAMTOOLS_FAIDX.out.fai, CREATE_REPEATS_BED.out.repeats_bed )
    ch_versions = ch_versions.mix(CREATE_REPMA_BED.out.versions)

    // Create the CpG bed file
    CREATE_CPG_BED ( reference )
    ch_versions = ch_versions.mix(CREATE_CPG_BED.out.versions)

    // Merge the CpG and repeat bed files
    CREATE_REPMA_CPG_BED (
        CREATE_REPMA_BED.out.genomefile,
        CREATE_REPMA_BED.out.sorted_repeats_bed,
        CREATE_REPMA_BED.out.ref_bed,
        CREATE_CPG_BED.out.cpg_bed)
    ch_versions = ch_versions.mix(CREATE_REPMA_CPG_BED.out.versions)

    emit:
    sorted_repeats_bed      = CREATE_REPMA_BED.out.sorted_repeats_bed           // channel: path(sorted_repeats_bed)
    repma_bed               = CREATE_REPMA_BED.out.repma_bed                    // channel: path(repma_bed)
    cpg_bed                 = CREATE_CPG_BED.out.cpg_bed                        // channel: path(cpg_bed)
    no_cpg_bed              = CREATE_REPMA_CPG_BED.out.no_cpg_bed               // channel: path(noCpG_bed)
    cpg_repeats_sorted_bed  = CREATE_REPMA_CPG_BED.out.cpg_repeats_sorted_bed   // channel: path(cpg_repeats_sorted_bed)
    no_cpg_repma_bed        = CREATE_REPMA_CPG_BED.out.no_cpg_repma_bed         // channel: path(noCpG_repma_bed)
    versions                = ch_versions                                       // channel: [ versions.yml ]
}