include { METASPADES; MEGAHIT; HIFIASM_META; METAFLYE } from '../modules/assembly'
include { CIRCULAR_CONTIGS_METAFLYE; CIRCULAR_CONTIGS_HIFIASM } from '../modules/circular_contigs.nf'
include { RENAME_CONTIGS } from '../modules/rename_contigs.nf'
include { QUAST as ASSEMBLY_QUAST} from '../modules/metaquast'
include { READS_DEDUPLICATION } from '../modules/reads_deduplication'
include { MINIMAP2 } from '../modules/read_alignment'
include { GET_ALIGNMENT_METRICS } from '../modules/read_alignment_manipulation'
include { MERGE_FASTQ } from '../modules/merge_fastq.nf'

workflow STEP_02_ASSEMBLY {
  take: 
  reads
  assembly
  has_assembly
  assembly_tool
  has_flowcell

  main:
    ch_assembler_v = Channel.empty() 
    ch_quast_v = Channel.empty()
    ch_bwa_mem_v = Channel.empty()
    ch_minimap2_v = Channel.empty()
    ch_samtools_v = Channel.empty()
    ch_infos_metaflye = Channel.empty()
    ch_circular = Channel.empty() 

    ch_assembly = assembly 
    ch_reads = reads

    if (!has_assembly & has_flowcell ){ 
			//////////////////
			// Manage Flowcell 
			//////////////////
			ch_reads_flowcell = reads
				.map { meta, fastq ->
					   [ meta.sample, meta, fastq ] }
				.groupTuple(by: [0])
				.branch {	id, meta, fastq ->
          single : fastq.size() == 1
            return [[id:meta.sample.unique().join(),
            sample:meta.sample.unique().join(),
            flowcell:meta.flowcell.join("_"),
            group:meta.group.unique().join(),
            assembly:meta.assembly.unique().join(),
            type:meta.type.unique().join()], fastq.flatten().sort{ it.baseName }  ]
          multiple: fastq.size() > 1
            return [[id:meta.sample.unique().join(),
            sample:meta.sample.unique().join(),
            flowcell:meta.flowcell.join("_"),
            group:meta.group.unique().join(),
            assembly:meta.assembly.unique().join(),
            type:meta.type.unique().join()], fastq.flatten().sort{ it.baseName }  ]
				}

			MERGE_FASTQ (ch_reads_flowcell.multiple)
				.reads
				.mix(ch_reads_flowcell.single)
				.set{ch_reads}
    }

    if (params.coassembly){ 
      ch_reads.map { meta, fastq ->
                   [ meta.group, meta, fastq] }
        .groupTuple(by: [0])
        .map { group, metas, fastq ->
              def meta = [:]
              meta.id = metas.group.unique().join()
              meta.sample = metas.sample.join("_")
              meta.flowcell = metas.flowcell.unique().join()
              meta.group =  metas.group.unique().join()
              meta.assembly = metas.assembly.unique().join()
              meta.type = metas.type.unique().join()
              if (params.type.toUpperCase() == "SR") {
                return [meta, fastq.collect { it[0] }, fastq.collect { it[1] }]
              } else {
                return [meta, fastq.flatten().sort{ it.baseName }]
              }}
        .set { ch_reads_assembly }

      if (has_assembly){ 
        ch_assembly = assembly.map { meta, assembly ->
                                   [ meta.group, meta, assembly] }
                              .groupTuple(by: [0]) 
                              .map{ group, metas, assembly ->
                                    def meta = [:]
                                    meta.id = metas.group.unique().join()
                                    meta.sample = metas.sample.join("_")
                                    meta.flowcell = metas.flowcell.unique().join()
                                    meta.group =  metas.group.unique().join()
                                    meta.assembly = metas.assembly.unique().join()
                                    meta.type = metas.type.unique().join()
                                    return [meta, assembly[0]] }
      }

    } else if (params.type.toUpperCase() == "SR") {
        ch_reads_assembly = ch_reads
        .map { meta, fastq ->
            return [meta, fastq[0], fastq[1]] }

    } else {
      ch_reads_assembly = ch_reads
    }

    if (!has_assembly){
      if (assembly_tool == 'metaspades') {
          METASPADES(ch_reads_assembly)
          ch_assembly = METASPADES.out.assembly
          ch_assembler_v = METASPADES.out.v_metaspades
      } else if (assembly_tool == 'megahit') {
          MEGAHIT(ch_reads_assembly)
          ch_assembly = MEGAHIT.out.assembly
          ch_assembler_v = MEGAHIT.out.v_megahit
      } else if (assembly_tool == 'hifiasm-meta') {
          HIFIASM_META(ch_reads_assembly)
          ch_assembly = HIFIASM_META.out.assembly
          ch_assembler_v = HIFIASM_META.out.v_hifiasm_meta
      } else if (assembly_tool == 'metaflye') {
          METAFLYE(ch_reads_assembly)
          ch_assembly = METAFLYE.out.assembly
          ch_infos_metaflye = METAFLYE.out.infos
          ch_assembler_v = METAFLYE.out.v_metaflye
      } else {
        exit 1,  "Invalid assembly parameter: ${assembly_tool}"
      }
    }

    RENAME_CONTIGS(ch_assembly)
    ch_assembly_renamed = RENAME_CONTIGS.out.fna

    ch_assembly_quast = ch_assembly_renamed
                        .map { meta, file -> file }
                        .collect(sort:{it.baseName})

    if (assembly_tool == 'hifiasm-meta') {
      CIRCULAR_CONTIGS_HIFIASM(ch_assembly_renamed)
      ch_circular = CIRCULAR_CONTIGS_HIFIASM.out.circular
    } else if (assembly_tool == 'metaflye') {
      ch_assembly_metaflye_infos = ch_assembly_renamed.join(ch_infos_metaflye)
      CIRCULAR_CONTIGS_METAFLYE(ch_assembly_metaflye_infos)
      ch_circular = CIRCULAR_CONTIGS_METAFLYE.out.circular
    }

    ASSEMBLY_QUAST( ch_assembly_quast,"02_assembly/02_1_primary_assembly/")
    ch_assembly_report = ASSEMBLY_QUAST.out.report
    ch_quast_v = ASSEMBLY_QUAST.out.v_quast

    ch_idxstats = Channel.empty()
    ch_flagstat = Channel.empty()
                                          
    if (params.coassembly){
      ch_reads.map { meta, fastq ->
                   [ meta.group, meta, fastq ]}
              .set { ch_reads_tmp }

      ch_assembly_renamed.map { meta, assembly ->
                              [ meta.group, assembly ]} 
                         .combine(ch_reads_tmp, by: 0)
                         .map { group, assembly, meta, fastq -> 
                              [ meta, assembly,fastq ]}
                         .set { ch_contigs_and_reads }
    } else {
      ch_contigs_and_reads = ch_assembly_renamed.join(ch_reads, remainder: true)
    }

    if ( params.type.toUpperCase() == "SR" ) {

      READS_DEDUPLICATION ( ch_contigs_and_reads )
      
      ch_reads = READS_DEDUPLICATION.out.dedup_reads
      ch_bam = READS_DEDUPLICATION.out.bam

      ch_bwa_mem_v = READS_DEDUPLICATION.out.v_bwa_mem2
      
    } else {
      MINIMAP2(ch_contigs_and_reads,"02_assembly/02_3_reads_vs_primary_assembly")
      ch_bam = MINIMAP2.out.bam

      ch_minimap2_v = MINIMAP2.out.v_minimap2
        
    }

    GET_ALIGNMENT_METRICS(ch_bam, "02_assembly/02_3_reads_vs_primary_assembly")

    ch_idxstats  = GET_ALIGNMENT_METRICS.out.sam_idxstat
    ch_flagstat =  GET_ALIGNMENT_METRICS.out.sam_flagstat
    ch_coverage =  GET_ALIGNMENT_METRICS.out.sam_coverage
    ch_samtools_v = GET_ALIGNMENT_METRICS.out.v_samtools

    ch_software_versions = ch_assembler_v.first().mix( ch_quast_v,
                                                       ch_bwa_mem_v.first(),
                                                       ch_minimap2_v.first(),
                                                       ch_samtools_v.first())


  emit:
    assembly = ch_assembly_renamed
    circular = ch_circular
    reads = ch_reads
    bam = ch_bam

    idxstats = ch_idxstats
    flagstat = ch_flagstat
    coverage = ch_coverage
    assembly_report = ch_assembly_report
    software_versions = ch_software_versions
}