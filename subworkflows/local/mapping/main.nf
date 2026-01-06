#! /usr/bin/env nextflow

// Indexing
include { SAMTOOLS_FAIDX as M_SAMTOOLS_FAIDX              } from '../../../modules/local/samtools/samtools_faidx.nf'
include { BWA_INDEX as M_BWA_INDEX                        } from '../../../modules/local/bwa/bwa_index.nf'
include { BOWTIE2_BUILD as M_BOWTIE2_BUILD                } from '../../../modules/local/bowtie2/bowtie2_build.nf'
include { BWA_ALN as M_BWA_ALN                            } from '../../../modules/local/bwa/bwa_aln.nf'
include { BWA_MEM as M_BWA_MEM                            } from '../../../modules/local/bwa/bwa_mem.nf'
include { SPLIT_FASTQ as M_SPLIT_FASTQ                    } from '../../../modules/local/awk/split_fastq.nf'
include { BWA_ALN as M_BWA_ALN_SHORT                      } from '../../../modules/local/bwa/bwa_aln.nf'
include { BWA_MEM as M_BWA_MEM_LONG                       } from '../../../modules/local/bwa/bwa_mem.nf'
include { SAMTOOLS_MERGE as M_BWA_ALN_MEM_MERGE           } from '../../../modules/local/samtools/samtools_merge.nf'
include { SAMTOOLS_MERGE as M_MERGED_UNMERGED_READS_BAM   } from '../../../modules/local/samtools/samtools_merge.nf'
include { BOWTIE2 as M_BOWTIE2                            } from '../../../modules/local/bowtie2/bowtie2.nf'
include { SAMTOOLS_INDEX as M_RAW_BAM_INDEX               } from '../../../modules/nf-core/samtools/index/main'


workflow MAPPING {
    take:
    reference
    reads

    main:
    ch_versions         = Channel.empty()

    M_SAMTOOLS_FAIDX ( reference )
    ch_versions         = ch_versions.mix(M_SAMTOOLS_FAIDX.out.versions)

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Index the reference genome if it is not already indexed
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Build the BWA index - only if BWA is selected as mapping tool
    def bwa_tools = ['bwa-aln', 'bwa-mem', 'bwa-aln-mem']
    if (bwa_tools.contains( params.mapping_tool_ancient ) || bwa_tools.contains( params.mapping_tool_modern )) {
        M_BWA_INDEX (reference, file(params.reference).getParent())
        ch_versions         = ch_versions.mix(M_BWA_INDEX.out.versions)
        ch_bwa_index        = M_BWA_INDEX.out.index_dir
    }

    // Build the Bowtie2 index - only if Bowtie2 is selected as mapping tool
    def bowtie2_tools = ['bowtie2']
    if (bowtie2_tools.contains( params.mapping_tool_ancient ) || bowtie2_tools.contains( params.mapping_tool_modern )) {
        M_BOWTIE2_BUILD (reference, file(params.reference).getParent())
        ch_versions         = ch_versions.mix(M_BOWTIE2_BUILD.out.versions)
        ch_bowtie2_index    = M_BOWTIE2_BUILD.out.index_dir
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. Mapping merged reads or single-end reads
    ////////////////////////////////////////////////////////////////////////////////////////////////

    ch_raw_bam = Channel.empty()

    // Branch reads into ancient and modern based on sample_type metadata
    reads.branch {
        ancient: it[0].sample_type == 'ancient'
        modern : it[0].sample_type == 'modern'
    }.set { ch_reads_branched }

    // Route reads to specific tools
    ch_reads_bwa_aln        = Channel.empty()
    ch_reads_bwa_mem        = Channel.empty()
    ch_reads_bwa_aln_mem    = Channel.empty()
    ch_reads_bowtie2        = Channel.empty()

    // Configure routing for ancient samples
    if (params.mapping_tool_ancient == 'bwa-aln')           { ch_reads_bwa_aln = ch_reads_bwa_aln.mix(ch_reads_branched.ancient) }
    else if (params.mapping_tool_ancient == 'bwa-mem')      { ch_reads_bwa_mem = ch_reads_bwa_mem.mix(ch_reads_branched.ancient) }
    else if (params.mapping_tool_ancient == 'bwa-aln-mem')  { ch_reads_bwa_aln_mem = ch_reads_bwa_aln_mem.mix(ch_reads_branched.ancient) }
    else if (params.mapping_tool_ancient == 'bowtie2')      { ch_reads_bowtie2 = ch_reads_bowtie2.mix(ch_reads_branched.ancient) }

    // Configure routing for modern samples
    if (params.mapping_tool_modern == 'bwa-aln')            { ch_reads_bwa_aln = ch_reads_bwa_aln.mix(ch_reads_branched.modern) }
    else if (params.mapping_tool_modern == 'bwa-mem')       { ch_reads_bwa_mem = ch_reads_bwa_mem.mix(ch_reads_branched.modern) }
    else if (params.mapping_tool_modern == 'bwa-aln-mem')   { ch_reads_bwa_aln_mem = ch_reads_bwa_aln_mem.mix(ch_reads_branched.modern) }
    else if (params.mapping_tool_modern == 'bowtie2')       { ch_reads_bowtie2 = ch_reads_bowtie2.mix(ch_reads_branched.modern) }


    // Run Mapping Tools
    // BWA ALN
    if (params.mapping_tool_ancient == 'bwa-aln' || params.mapping_tool_modern == 'bwa-aln') {
        M_BWA_ALN ( ch_reads_bwa_aln, ch_bwa_index )
        ch_versions         = ch_versions.mix(M_BWA_ALN.out.versions)
        ch_raw_bam          = ch_raw_bam.mix(M_BWA_ALN.out.bam)
    }
    // BWA MEM
    if (params.mapping_tool_ancient == 'bwa-mem' || params.mapping_tool_modern == 'bwa-mem') {
        M_BWA_MEM ( ch_reads_bwa_mem, ch_bwa_index )
        ch_versions         = ch_versions.mix(M_BWA_MEM.out.versions)
        ch_raw_bam          = ch_raw_bam.mix(M_BWA_MEM.out.bam)
    }
    // BWA ALN-MEM
    if (params.mapping_tool_ancient == 'bwa-aln-mem' || params.mapping_tool_modern == 'bwa-aln-mem') {
        //split fastq
        M_SPLIT_FASTQ ( ch_reads_bwa_aln_mem )
        ch_versions         = ch_versions.mix(M_SPLIT_FASTQ.out.versions)
        //align short reads with BWA ALN
        M_BWA_ALN_SHORT ( M_SPLIT_FASTQ.out.short_reads, ch_bwa_index )
        ch_versions         = ch_versions.mix(M_BWA_ALN_SHORT.out.versions)
        //align long reads with BWA MEM
        M_BWA_MEM_LONG ( M_SPLIT_FASTQ.out.long_reads,  ch_bwa_index )
        ch_versions         = ch_versions.mix(M_BWA_MEM_LONG.out.versions)
        //merge BAMs from short and long reads
        ch_raw_bam_aln_mem  = M_BWA_ALN_SHORT.out.bam.join(M_BWA_MEM_LONG.out.bam)
            .map { meta, file1, file2 -> [meta, [file1, file2]] }
        M_BWA_ALN_MEM_MERGE ( ch_raw_bam_aln_mem, reference )
        ch_versions         = ch_versions.mix(M_BWA_ALN_MEM_MERGE.out.versions)
        ch_raw_bam          = ch_raw_bam.mix(M_BWA_ALN_MEM_MERGE.out.bam)
    }
    // BOWTIE2
    if (params.mapping_tool_ancient == 'bowtie2' || params.mapping_tool_modern == 'bowtie2') {
        M_BOWTIE2 ( ch_reads_bowtie2, ch_bowtie2_index )
        ch_versions         = ch_versions.mix(M_BOWTIE2.out.versions)
        ch_raw_bam          = ch_raw_bam.mix(M_BOWTIE2.out.bam)
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 3. Merge the mapped unmerged reads if provided
    ////////////////////////////////////////////////////////////////////////////////////////////////

    def ch_merged_raw_bam = null
    if (params.merge_reads.toBoolean() && params.keep_unmerged_reads.toBoolean()) {
        // Group the unmerged reads BAMs by sample ID (removing the '-unmerged' suffix)
        ch_raw_bam_grouped = ch_raw_bam.map { meta, bam ->
                def new_meta = meta.clone()
                new_meta.id = meta.id.replace("-unmerged", "")
                tuple(new_meta.id, new_meta, bam)
            }
            .groupTuple()
            .map { id, metas, bams ->
                // pick one meta; choose the one with single_end == false if present
                def final_meta = metas.find { !it.single_end } ?: metas[0]
                // remove grouping id, return meta + joined bam list
                tuple(final_meta, bams)
            }

        M_MERGED_UNMERGED_READS_BAM ( ch_raw_bam_grouped, reference )
        ch_versions = ch_versions.mix(M_MERGED_UNMERGED_READS_BAM.out.versions)
        ch_merged_raw_bam = M_MERGED_UNMERGED_READS_BAM.out.bam
    }

    // Use merged raw bam if created, else use original raw bam
    ch_raw_bam_final = ch_merged_raw_bam ?: ch_raw_bam

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 4. Index the raw BAM files
    ////////////////////////////////////////////////////////////////////////////////////////////////

    M_RAW_BAM_INDEX ( ch_raw_bam_final )
    ch_versions = ch_versions.mix(M_RAW_BAM_INDEX.out.versions)

    ////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    fai         = M_SAMTOOLS_FAIDX.out.fai       // channel: path(index)
    raw_bam     = ch_raw_bam_final               // channel: [ val(meta), [ bam ] ]
    raw_bai     = M_RAW_BAM_INDEX.out.bai        // channel: [ val(meta), [ bai ] ]
    versions    = ch_versions                    // channel: [ versions.yml ]
}