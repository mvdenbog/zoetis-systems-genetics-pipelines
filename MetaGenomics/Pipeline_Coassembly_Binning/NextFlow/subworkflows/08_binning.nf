include { GENERATE_DEPTH_FILES; METABAT2; MAXBIN2; CONCOCT; BINETTE; UNBINNED_CONTIGS } from '../modules/binning'
include { CHECKM2 } from '../modules/checkm2'
include { DREP } from '../modules/drep'
include { BWA_MEM;BWA_MEM as BWA_MEM_CROSS_ALIGNMENT; MINIMAP2; MINIMAP2 as MINIMAP2_CROSS_ALIGNMENT } from '../modules/read_alignment'
include { GET_ALIGNMENT_METRICS } from '../modules/read_alignment_manipulation'
include { GTDBTK } from '../modules/gtdbtk'
include { GENOMES_ABUNDANCES_PER_SAMPLE; ADD_QUAST_INFO_TO_BINS; BINS_STATS_TO_MUTLIQC_FORMAT} from '../modules/sum_up_bins_informations'

// Declare syntax version
nextflow.enable.dsl=2
// Script parameters

workflow STEP_08_BINNING {
  take: reads
  assembly
  bam
  gtdbtk_db
  checkm2_db
  quast
  circular

  main:

  ch_bwa_mem_v = Channel.empty()
  ch_minimap2_v = Channel.empty()
  ch_samtools_v = Channel.empty()
  ch_metabat2_v = Channel.empty()
  ch_maxbin_v = Channel.empty()
  ch_concoct_v = Channel.empty()
  ch_binette_v = Channel.empty()
  ch_checkm2_v = Channel.empty()
  ch_dRep_v = Channel.empty()
  ch_gtdbtk_v = Channel.empty()
  ch_bins_abundances_report = Channel.empty()
  ch_bins_stats_report = Channel.empty()
  ch_metabat_bins = Channel.empty()
  ch_maxbin_bins = Channel.empty()
  ch_concoct_bins = Channel.empty()


  ////// 
  /// Config
  //////
  ch_heatmap_header_multiqc = file("$projectDir/assets/multiqc/bins_abundances_header.txt", checkIfExists: true)
  ch_barplot_bin_count_header_multiqc = file("$projectDir/assets/multiqc/barplot_bin_count_header.txt", checkIfExists: true)
  ch_barplot_bin_size_header_multiqc = file("$projectDir/assets/multiqc/barplot_bin_size_header.txt", checkIfExists: true)
  ch_table_header_multiqc = file("$projectDir/assets/multiqc/bin_stats_table_header.txt", checkIfExists: true)
  
  ////// 
  /// CROSS ALIGNMENT 
  //////

  if (params.binning_cross_alignment != 'individual'){
    if (params.binning_cross_alignment == 'all') {
      // combine assemblies with reads of all samples
      ch_reads_assembly = reads.combine(assembly)
                          .map{ meta_reads, reads, meta_assembly, assembly ->
                            if (((meta_reads.id != meta_assembly.id) && !(params.coassembly)) || (params.coassembly && (meta_reads.group != meta_assembly.group))){
                              [[id:meta_reads.id+"_"+meta_assembly.id, sample:meta_assembly.sample, flowcell:meta_assembly.flowcell, group:meta_assembly.group, assembly:meta_assembly.assembly, type:meta_assembly.type], assembly, reads]
                            }
                          }
    } else if (params.binning_cross_alignment == 'group'){
      // combine assemblies with reads of samples from same group
      ch_assembly_group = assembly.map{ meta, assembly -> [ meta.group, meta, assembly ] }
      ch_reads_assembly = reads.map{ meta, reads -> [ meta.group, meta, reads ] }
                          .combine(ch_assembly_group, by: 0)
                          .map{ group, meta_reads, reads, meta_assembly, assembly ->
                            if (meta_reads != meta_assembly){
                              [[id:meta_reads.id+"_"+meta_assembly.id, sample:meta_assembly.sample, flowcell:meta_assembly.flowcell, group:meta_assembly.group, assembly:meta_assembly.assembly, type:meta_assembly.type], assembly, reads]
                            }
                          }
    }
    // cross alignment
    if (params.type == 'SR') {
      BWA_MEM_CROSS_ALIGNMENT(ch_reads_assembly,"08_binning/08_1_binning_per_sample/cross_mapping")
      ch_bam_cross_alignment = BWA_MEM_CROSS_ALIGNMENT.out.bam
    } else {
      MINIMAP2_CROSS_ALIGNMENT(ch_reads_assembly,"08_binning/08_1_binning_per_sample/cross_mapping")
      ch_bam_cross_alignment = MINIMAP2_CROSS_ALIGNMENT.out.bam
    } 
    // formatting channel  
    ch_bam = ch_bam_cross_alignment.mix(bam)
            .map { meta, bam, bai -> 
                  if (params.coassembly){[ meta.group, meta, bam, bai ]} 
                  else  {[ meta.sample, meta, bam, bai ]}}
            .groupTuple(by: [0])
            .map { sample,metas, bam, bai ->
              [ metas.min { it.id.size() }, bam, bai ]
            } 
  } else {
      ch_bam = bam
  }

  if (params.coassembly){
    if (params.binning_cross_alignment == 'all') {
      ch_bam.map { meta, bam, bai -> [ meta.group, bam, bai ]  }
            .set { ch_bam_tmp }
    } else {
      ch_bam.map { meta, bam, bai -> [ meta.group, meta, bam, bai ]  }
            .groupTuple(by: [0])
            .map { group, metas, bam, bai ->
                [ group, bam, bai ]} 
            .set { ch_bam_tmp }
    }
    assembly.map { meta, assembly ->
                  [ meta.group, meta, assembly ]} 
            .combine( ch_bam_tmp, by: 0)
            .map { group, meta, assembly, bam, bai -> 
                  [ meta, assembly, bam, bai]}
            .tap { ch_assembly_bam }
            .map { meta, assembly, bam, bai ->
                  [ meta, assembly ]}
            .tap { ch_assembly }

    ch_assembly_bam.map{  meta, assembly, bam, bai ->
                       [ meta, bam, bai ]}
                   .set{ ch_bam}
  } else {
    ch_assembly_bam = assembly.join(ch_bam)
    ch_assembly = assembly
  }
  ///////////
  /// BINNING
  ///////////

  ch_depth = GENERATE_DEPTH_FILES(ch_bam)
  ch_assembly_depth = ch_assembly.join(ch_depth)

  METABAT2(ch_assembly_depth)
  ch_metabat_bins = METABAT2.out.bins.filter{ t -> t[1].list().size()}
  ch_metabat2_v = METABAT2.out.v_metabat2

  MAXBIN2(ch_assembly_depth)
  ch_maxbin_bins = MAXBIN2.out.bins.filter{ t -> t[1].list().size()}
  ch_maxbin_v = MAXBIN2.out.v_maxbin

  CONCOCT(ch_assembly_bam)
  ch_concoct_bins = CONCOCT.out.bins.filter{ t -> t[1].list().size()}
  ch_concoct_v = CONCOCT.out.v_concoct

  ch_circular = circular.filter{ t -> t[1].list().size()}
  //////////////////
  //// BIN REFINEMENT  
  //////////////////
  ch_bins_set = ch_metabat_bins.mix(ch_concoct_bins, ch_maxbin_bins, ch_circular)
                                .groupTuple(by: [0])
                                .branch{ meta, bins ->
                                    multiple: bins.size() > 1
                                    single: bins.size() == 1
                                }

  CHECKM2(ch_bins_set.single, checkm2_db)
  ch_checkm2_v = CHECKM2.out.v_checkm2


  ch_all_bins_stats = Channel.empty()
  ch_binette_stats_and_quast = Channel.empty()
  ch_binette_bins = Channel.empty()
  ch_binette_stats = Channel.empty()
  ch_bins_assembly = Channel.empty()
  ch_drep_fna = Channel.empty()
  ch_binette_bins_filter = Channel.empty()
  ch_binette_stats_filter = Channel.empty()
  ch_bins_drep = Channel.empty()
  ch_bam_bins = Channel.empty()
  ch_reads_fna = Channel.empty()
  ch_gtdbtk_affi = Channel.empty()
  ch_drep_stats = Channel.empty()

  ch_bins_assembly = ch_bins_set.multiple.join(ch_assembly)

  BINETTE(ch_bins_assembly, params.min_completeness, params.max_contamination, checkm2_db)
  
  ch_binette_bins = BINETTE.out.bins

  ch_binette_stats = BINETTE.out.checkm_stats.join(ch_binette_bins)
                                              .map { meta, stats, bins -> 
                                                    [ meta, stats ] }
  
  ch_bins_assembly = ch_binette_bins.join(ch_assembly)

  UNBINNED_CONTIGS(ch_bins_assembly)

  ch_binette_stats_and_quast =  ch_binette_stats.combine(quast) 
  ch_binette_v = BINETTE.out.v_binette
  
  ADD_QUAST_INFO_TO_BINS(ch_binette_stats_and_quast)
  
  ch_all_bins_stats = ADD_QUAST_INFO_TO_BINS.out.bins_stats.collect()

  BINS_STATS_TO_MUTLIQC_FORMAT(ch_all_bins_stats, 
                                  ch_barplot_bin_count_header_multiqc, 
                                  ch_barplot_bin_size_header_multiqc)
  ch_bins_stats_report = BINS_STATS_TO_MUTLIQC_FORMAT.out.report
  
  ////////////////////////
  ///// DEREPLICATION BIN 
  ////////////////////////

  ch_bins_filter=ch_binette_bins.filter{ it[1].getClass() == java.util.ArrayList }
  ch_binette_bins_filter = ch_bins_filter.map { meta, file -> file}
                                         .collect()

  ch_binette_stats_filter = ch_binette_stats.join(ch_bins_filter)
                                            .map { meta, stats, bins -> stats }
                                            .collect()


  DREP(ch_binette_stats_filter, ch_binette_bins_filter, params.drep_threshold)
  ch_bins_drep = DREP.out.drep_bins_folder
  ch_dRep_v = DREP.out.v_dRep
  ch_drep_fna = DREP.out.fna
  ch_drep_stats = DREP.out.output_drep_stats

  ///////////////////////////////
  ///// TAXONOMIC AFFILIATION BIN
  /////////////////////////////// 

  GTDBTK(ch_bins_drep, gtdbtk_db)
  ch_gtdbtk_v = GTDBTK.out.v_gtdbtk
  ch_gtdbtk_affi = GTDBTK.out.gtdbtk_affiliations_predictions

  /////////////////////////////
  ////GENOMES ABUNDANCES
  /////////////////////////////
  ch_reads_fna = reads.combine(ch_drep_fna)
                .map { meta, reads, bins -> [meta, bins, reads] }

  if (params.type == 'SR') {
    BWA_MEM(ch_reads_fna, "08_binning/08_4_mapping_on_final_bins/mapping")
    ch_bam_bins = BWA_MEM.out.bam
    ch_bwa_mem_v = BWA_MEM.out.v_bwa_mem2
    ch_samtools_v = BWA_MEM.out.v_samtools
  } else {
    MINIMAP2(ch_reads_fna, "08_binning/08_4_mapping_on_final_bins/mapping")
    ch_bam_bins = MINIMAP2.out.bam
    ch_minimap2_v = MINIMAP2.out.v_minimap2
    ch_samtools_v = MINIMAP2.out.v_samtools
  }

  ch_collect_coverages = Channel.empty()
  ch_collect_flagstats = Channel.empty()
  GET_ALIGNMENT_METRICS(ch_bam_bins, "08_binning/08_4_mapping_on_final_bins/stats")

  ch_collect_coverages = GET_ALIGNMENT_METRICS.out.sam_coverage.map {id, file -> file}
              .collect()
  ch_collect_flagstats = GET_ALIGNMENT_METRICS.out.sam_flagstat.collect()

  GENOMES_ABUNDANCES_PER_SAMPLE(ch_collect_coverages, ch_collect_flagstats, 
   ch_bins_drep, ch_drep_stats , ch_gtdbtk_affi ,
   ch_heatmap_header_multiqc, ch_table_header_multiqc)

  ch_bins_abundances_report = GENOMES_ABUNDANCES_PER_SAMPLE.out.report

  ch_software_versions = ch_bwa_mem_v.first().mix(ch_minimap2_v.first(),
                                                  ch_samtools_v.first(),
                                                  ch_metabat2_v.first(),
                                                  ch_maxbin_v.first(),
                                                  ch_concoct_v.first(),
                                                  ch_binette_v,
                                                  ch_checkm2_v.first(),
                                                  ch_dRep_v,
                                                  ch_gtdbtk_v
                                                  )

  emit:
    bins_abundances_report = ch_bins_abundances_report
    bins_stats_report = ch_bins_stats_report
    software_versions = ch_software_versions

}
