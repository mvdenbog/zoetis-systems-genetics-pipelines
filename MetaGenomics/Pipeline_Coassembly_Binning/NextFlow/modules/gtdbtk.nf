process GTDBTK {
 publishDir "${params.outdir}/08_binning/08_3_gtdbtk", mode: 'copy'
 label 'BINNING'
  input:
    val(drep_bins_folder)
    val(gtdbtk_db)
      
  output:
    path("gtdbtk.bac120.summary.tsv*"), emit : gtdbtk_affiliations_predictions
    path("v_gtdbtk.txt"),emit : v_gtdbtk
  script:
  """

  export GTDBTK_DATA_PATH=$gtdbtk_db

  gtdbtk classify_wf --genome_dir $drep_bins_folder -x fa --out_dir ./ --pplacer_cpus ${task.cpus} --cpus ${task.cpus}
  echo \$(gtdbtk -h 2>&1) &> v_gtdbtk.txt
  """
}