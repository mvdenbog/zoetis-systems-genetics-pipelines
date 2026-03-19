process DREP {
 publishDir "${params.outdir}/08_binning/08_2_dereplicated_bins", mode: 'copy', pattern: 'dereplicated_ge*'
 publishDir "${params.outdir}/08_binning/08_2_dereplicated_bins", mode: 'copy', pattern: 'figur*'
 publishDir "${params.outdir}/08_binning/08_2_dereplicated_bins", mode: 'copy', pattern: 'data_tables/*'
 label 'BINNING'
  input:
    path "stats"
    path "drep_input_bins/*"
    val drep_threshold
      
  output:
    path("dereplicated_genomes/"), emit : drep_bins_folder
    path("figures/"), emit : figures
    path("data_tables/genomeInformation.csv"), emit : output_drep_stats
    path("data_tables/Bins_clusters_composition.tsv"), emit : bins_clusters_composition
    path("dereplicated_bins.fna"), emit : fna
    path("v_dRep.txt"), emit : v_dRep

  script:
  """

  echo "genome,completeness,contamination" > ./bins_stats.csv
  cut -f1,2,3 $stats | awk -F "\t" '{if (\$1 ~ /_bin./ ) {print \$1".fa,"\$2","\$3}}' >> ./bins_stats.csv

  dRep dereplicate ./ -p ${task.cpus} --S_algorithm fastANI -sa ${drep_threshold} -g drep_input_bins/* --genomeInfo bins_stats.csv \
  --completeness 0 --contamination 100 --length 5000

  retrieve_bins_cluster_composition.py --clusters ./data_tables/Cdb.csv --representatives ./data_tables/Wdb.csv -o ./data_tables/Bins_clusters_composition.tsv

  for fi in \$(ls ./dereplicated_genomes/*)
  do
    cat \$fi >> ./dereplicated_bins.fna
  done

  echo \$(dRep -h 2>&1) > v_dRep.txt
  """
}