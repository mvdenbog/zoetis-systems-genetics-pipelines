include { CHUNK_ASSEMBLY_FILTER; MERGE_ASSEMBLY_FILTER} from '../modules/filtering_cpm.nf'
include { QUAST } from '../modules/metaquast'
include { MINIMAP2; BWA_MEM } from '../modules/read_alignment'
include { GET_ALIGNMENT_METRICS} from '../modules/read_alignment_manipulation'

workflow STEP_03_ASSEMBLY_FILTER {
  take:
    assemblies
    reads
    idxstats
    bam
    min_cpm
    

  main:

    filtering_outdir = "03_filtering/filtering_at_${params.min_contigs_cpm}cpm"

    // if filtering has no effect on assembly. A symblink will be created between reads vs assembly of the step 02 in the step 03 outdir
    unfiltered_assembly_bam_outdir = "02_assembly/02_3_reads_vs_primary_assembly" 
    filtered_assembly_bam_outdir = "${filtering_outdir}/03_2_reads_vs_filtered_assembly/"
    

    ch_chunk_assembly_for_filter = assemblies
                                  .splitFasta(by: 100000, file: true)
// .map { meta, fasta ->
//     tuple(meta, fasta)
// }
// .map { tuple ->
//     def meta  = tuple[0]
//    def fasta = tuple[1]
//    meta.__chunk_idx = (meta.__chunk_idx ?: 0) + 1
//     def idx = meta.__chunk_idx
//    [ meta, fasta, idx ]
// }
// .map { meta, fasta, idx ->
//     def name = "${meta.id}_chunk${idx}"
//    fasta = fasta.copyTo("${params.outdir}/chunks/${name}.fasta")
//    [ meta, fasta ]
//}


    if (params.coassembly){ 
      idxstats.map { meta, idxstats ->
                   [ meta.group, meta, idxstats] }
              .groupTuple(by: [0])
              .map { group, metas, idxstats -> 
                    def meta = [:]
                    meta.id = metas.group.unique().join()
                    meta.sample = metas.sample.join("_")
                    meta.flowcell = metas.flowcell.unique().join()
                    meta.group =  metas.group.unique().join()
                    meta.assembly = metas.assembly.unique().join()
                    meta.type = metas.type.unique().join()
                    [ group, meta, idxstats] }
              .set { ch_idxstats_tmp }
      ch_chunk_assembly_for_filter.map { meta, assembly ->
                                       [ meta.group, assembly ]} 
                                  .combine(ch_idxstats_tmp, by: 0)
                                  .map { group, assembly, meta, idxstats -> 
                                       [ meta, assembly, idxstats ]}
                                  .set { ch_assembly_and_idxstats }
    } else {
      ch_assembly_and_idxstats = ch_chunk_assembly_for_filter
                                .combine(idxstats, by:0)
    }

    

    CHUNK_ASSEMBLY_FILTER (
      ch_assembly_and_idxstats,
      min_cpm
    )
    ch_chunk_selected =  CHUNK_ASSEMBLY_FILTER.out.chunk_selected
    ch_chunk_discarded = CHUNK_ASSEMBLY_FILTER.out.chunk_discarded

    ch_chunk_selected
      .groupTuple(by: 0)
      .set{ ch_grouped_selected }

    ch_chunk_discarded
      .groupTuple(by: 0)
      .set{ ch_grouped_discarded }

    MERGE_ASSEMBLY_FILTER (
      ch_grouped_selected,
      ch_grouped_discarded,
      min_cpm,
      "${filtering_outdir}/03_1_filtered_assembly"
    )
    ch_merged_selected = MERGE_ASSEMBLY_FILTER.out.merged_selected
    discarded_contigs = MERGE_ASSEMBLY_FILTER.out.merged_discarded

    ch_merged_selected_all = ch_merged_selected
                            .map { meta, file -> file }
                            .collect(sort:{it.baseName})
    QUAST( ch_merged_selected_all, "${filtering_outdir}/03_1_filtered_assembly/" )
    ch_quast_report = QUAST.out.report


    // Differentiate sample with and without discarded_contigs 
    // samples with no discarded_contigs are not going to be processed to process 
    discarded_contigs_category = discarded_contigs.branch{ 
                                  without: it[1].isEmpty() 
                                  with: !it[1].isEmpty() 
                                  }

  
    if (params.coassembly){
      discarded_contigs_category.without.map { meta, discarded_empty -> [ meta.group ]} 
                                        .combine( bam.map { meta, bam, bai ->
                                                          [ meta.group, meta, bam, bai ]}, by: 0)
                                        .map{ group, meta, bam, bai -> 
                                            [ meta, bam, bai ]}
                                        .set{ ch_bam_unchanged_by_filtering }
                                      
      ch_selected_contigs_and_reads= discarded_contigs_category.with.map {meta, discarded_contigs -> meta.group}
                                                                    .join( ch_merged_selected.map { meta, contigs -> 
                                                                                                  [meta.group, meta, contigs]})
                                                                    .combine( reads.map{ meta, reads ->
                                                                                       [ meta.group, meta, reads ]}, by: 0)
                                                                    .map{ group, meta_contigs, contigs, meta_reads, reads -> 
                                                                        [ meta_reads, contigs, reads ]}

    } else {
      ch_bam_unchanged_by_filtering =  discarded_contigs_category.without.map{ it -> it[0]}
                                                                        .join(bam)

      ch_selected_contigs_and_reads = discarded_contigs_category.with.map{ it -> it[0]}
                                                                    .join(ch_merged_selected).join(reads)
    }

    // make a symblink with the bam and bai from step 02 for samples that have not been affected by the filtering (no contig discarded)
    result_path_dir = file("${params.outdir}/${filtered_assembly_bam_outdir}/")
    result_path_dir.mkdirs()
    
    ch_bam_unchanged_by_filtering.map { meta, bam, bai ->
         { file("${result_path_dir}/${meta.id}/").mkdir()
           file("${params.outdir}/${unfiltered_assembly_bam_outdir}/${meta.id}/${meta.id}.bam")
                .mklink("${result_path_dir}/${meta.id}/${meta.id}.bam", overwrite:true)
           file("${params.outdir}/${unfiltered_assembly_bam_outdir}/${meta.id}/${meta.id}.bam.bai")
                .mklink("${result_path_dir}/${meta.id}/${meta.id}.bam.bai", overwrite:true)
      }
    }

    if ( params.type.toUpperCase() == "SR" ) {
      BWA_MEM(ch_selected_contigs_and_reads, filtered_assembly_bam_outdir)
      ch_bam_post_filtering = BWA_MEM.out.bam
    }
    else {
      MINIMAP2(ch_selected_contigs_and_reads, filtered_assembly_bam_outdir)
      ch_bam_post_filtering = MINIMAP2.out.bam
    }


    ch_all_bam = ch_bam_unchanged_by_filtering.mix(ch_bam_post_filtering)

    GET_ALIGNMENT_METRICS(ch_all_bam, filtered_assembly_bam_outdir)


    ch_flagstat =  GET_ALIGNMENT_METRICS.out.sam_flagstat
    ch_coverage =  GET_ALIGNMENT_METRICS.out.sam_coverage

  


    emit:
      selected_contigs = ch_merged_selected
      quast_report = ch_quast_report
      bam = ch_all_bam
      sam_coverage = ch_coverage
      sam_flagstat = ch_flagstat

}

