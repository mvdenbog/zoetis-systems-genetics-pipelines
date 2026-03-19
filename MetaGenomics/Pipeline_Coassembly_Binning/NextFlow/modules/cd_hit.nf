process INDIVIDUAL_CD_HIT {
  tag "${meta.id}"
  publishDir "${params.outdir}/06_func_annot/06_1_clustering", mode: 'copy', pattern: "${meta.id}.cd-hit-est*"
  label 'CD_HIT'

  input:
  tuple val(meta), path(ffn)
  val pct_id

  output:
  path("${meta.id}.cd-hit-est.${pct_id}.fasta"), emit: clstr_fasta
  path("${meta.id}.cd-hit-est.${pct_id}.table_cluster_contigs.txt"), emit: clstr_table
  path "v_cdhit.txt", emit: v_cdhit

  script:
  """
  cd-hit-est -c ${pct_id} -i ${ffn} -o ${meta.id}.cd-hit-est.${pct_id}.fasta -T ${task.cpus} -M ${task.mem} -d 150
  cd_hit_produce_table_clstr.py -i ${meta.id}.cd-hit-est.${pct_id}.fasta.clstr -o ${meta.id}.cd-hit-est.${pct_id}.table_cluster_contigs.txt
  echo \$(cd-hit -h 2>&1) > v_cdhit.txt
  """
}



// Global clustering with CD-HIT.
process GLOBAL_CD_HIT {
  publishDir "${params.outdir}/06_func_annot/06_1_clustering", mode: 'copy'
  label 'CD_HIT'

  input:
  path cluster_fasta 
  val pct_id

  output:
  path "All-cd-hit-est.${pct_id}.fasta", emit: fasta_clusters
  path "table_clstr.txt", emit: clstr_table

  script:
  """
  # *fasta is important to get the correct order
  cat *.fasta > All-cd-hit-est.${pct_id}
  cd-hit-est -c ${pct_id} -i All-cd-hit-est.${pct_id} -o All-cd-hit-est.${pct_id}.fasta -T ${task.cpus} -M ${task.mem} -d 150
  cd_hit_produce_table_clstr.py -i All-cd-hit-est.${pct_id}.fasta.clstr -o table_clstr.txt
  """
}

workflow CD_HIT {
  
take:
ch_assembly // channel: [ val(meta.id), path(assemblyfasta) ]
ch_percentage_identity // channel: val

main:
  INDIVIDUAL_CD_HIT( ch_assembly, ch_percentage_identity )
  ch_individual_clusters = INDIVIDUAL_CD_HIT.out.clstr_fasta.collect()
  GLOBAL_CD_HIT(ch_individual_clusters , ch_percentage_identity )

emit:
  individual_clstr_table = INDIVIDUAL_CD_HIT.out.clstr_table
  global_clstr_table = GLOBAL_CD_HIT.out.clstr_table
  v_cdhit = INDIVIDUAL_CD_HIT.out.v_cdhit

}
