#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.skip_tools = 'markduplicates'
params.split_fastq = -1
params.save_reference = true
params.save_output_as_bam = false


include { BWAMEM2_INDEX                          } from './modules/bwamem2/index/main'
include { BWAMEM2_MEM            } from './modules/bwamem2/mem/main'
include { SAMTOOLS_FAIDX                         } from './modules/samtools/faidx/main'
include { SAMTOOLS_INDEX as INDEX_MERGE_BAM } from './modules/samtools/index/main'
include { FASTQC                                      } from './modules/fastqc/main'

// Convert BAM files
include { SAMTOOLS_CONVERT as BAM_TO_CRAM             } from './modules/samtools/convert/main'
// Merge and index BAM files (optional)
include { BAM_MERGE_INDEX_SAMTOOLS                    } from './subworkflows/bam_merge_index_samtools/main'
include { BAM_VARIANT_CALLING_DEEPVARIANT                                              } from './subworkflows/bam_variant_calling_deepvariant/main'
include { GATK4_CREATESEQUENCEDICTIONARY         } from './modules/gatk4/createsequencedictionary/main'
// Build intervals if needed
include { PREPARE_INTERVALS                           } from './subworkflows/prepare_intervals/main'
include { DEEPVARIANT                               } from './modules/deepvariant/main'
include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_GVCF } from './modules/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_VCF  } from './modules/gatk4/mergevcfs/main'

workflow {

params.intervals = "" 

Channel.fromFilePairs(params.input, flat: false)
.map { sample,  reads ->
        meta = [
            sample:sample,
        ]
        [ "sample":sample+"_"+reads.first().name.split('_')[1],
          "data_type":"fastq",
          "read_group":"\"@RG\tID:"+sample+"\tSM:"+sample+"\tLB:LB\tPU:PU\tPL:Illumina\"",
          "fastq_1":reads.first(),
          "fastq_2":reads.last()
        ]
    }
.set { ch_inputs}

// ch_inputs.view()

 ch_inputs
    .map { item ->
        def meta = [:]
        meta.id = item.sample
        meta.read_group = item.read_group
        meta.data_type = item.data_type
        meta.sample = item.sample
          return [meta,[item.fastq_1,item.fastq_2]]
      }
    .set{input_sample}

fasta_ch = Channel
     .fromPath(params.fasta)

fasta = fasta_ch.map{ fasta_ch -> [ [ id:fasta_ch.baseName ], fasta_ch ] }

SAMTOOLS_FAIDX(fasta, [['id':null], []])

                if ( params.host_index == "" ||!params.host_index ) {

BWAMEM2_INDEX(fasta_ch)
index_alignment = BWAMEM2_INDEX.out.index.map{ meta, index -> [index] }.collect()
                }
                else {
                    index_alignment = Channel.value(file(params.host_index))
                }



fasta_fai             = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }

input_fastq = input_sample

FASTQC(input_fastq)

sort_bam = true
BWAMEM2_MEM(input_fastq, index_alignment.map{ it -> [ [ id:'index' ], it ] }, sort_bam)
bam_mapped = BWAMEM2_MEM.out.bam

INDEX_MERGE_BAM(bam_mapped)

bam_bai = bam_mapped.join(INDEX_MERGE_BAM.out.bai)

GATK4_CREATESEQUENCEDICTIONARY(fasta)

dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict

PREPARE_INTERVALS(fasta_fai, params.intervals, params.no_intervals)
intervals= PREPARE_INTERVALS.out.intervals_bed  

intervals_and_num_intervals = intervals.map{ interval, num_intervals ->
        if ( num_intervals < 1 ) [ [], num_intervals ]
        else [ interval, num_intervals ]
    }

// Combine cram and intervals for spread and gather strategy
cram_intervals = bam_bai.combine(intervals)
// Move num_intervals to meta map
.map{ meta, cram, crai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram, crai, intervals ]}

DEEPVARIANT(cram_intervals, fasta_ch.map{ fasta_ch -> [ [ id:fasta_ch.baseName ], fasta_ch ] }.first(), fasta_fai.map{ fasta_fai -> [ [ id:fasta_fai.baseName ], fasta_fai ] }.first(), [ [ id:'null' ], [] ])

// Figuring out if there is one or more vcf(s) from the same sample
    vcf_out = DEEPVARIANT.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more gvcf(s) from the same sample
    gvcf_out = DEEPVARIANT.out.gvcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }


    // Only when using intervals
    gvcf_to_merge = gvcf_out.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    vcf_to_merge = vcf_out.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()

    MERGE_DEEPVARIANT_GVCF(gvcf_to_merge, dict.first())
    MERGE_DEEPVARIANT_VCF(vcf_to_merge, dict.first())

    // Mix intervals and no_intervals channels together
    gvcf = Channel.empty().mix(MERGE_DEEPVARIANT_GVCF.out.vcf, gvcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'deepvariant' ], vcf ] }

    // Mix intervals and no_intervals channels together
    vcf = Channel.empty().mix(MERGE_DEEPVARIANT_VCF.out.vcf, vcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'deepvariant' ], vcf ] }

}

