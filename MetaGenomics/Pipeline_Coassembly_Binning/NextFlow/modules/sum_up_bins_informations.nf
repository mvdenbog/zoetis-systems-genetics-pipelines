process GENOMES_ABUNDANCES_PER_SAMPLE {
   publishDir "${params.outdir}/08_binning", mode: 'copy'
   
   input:
      path coverage_files
      path flagstats_files
      val(bins_folder)
      path genomes_informations
      path affiliations_predictions
      path heatmap_header_mqc
      path table_header_mqc

   output:
      path "genomes_abundances.tsv" , emit: genomes_abundances
      tuple path("stats/genomes_abundances_mqc.tsv"), path("stats/genomes_checkm_mqc.json"), path("stats/bins_general_stats_mqc.tsv"), emit: report

   script:
      """
      mkdir -p stats
      bins_per_sample_summarize.py --list_of_coverage_files ${coverage_files} \
      --list_of_flagstats_files ${flagstats_files} --affiliations_predictions ${affiliations_predictions} \
      --bins_folder ${bins_folder} --genomes_informations ${genomes_informations} \
      --output_file genomes_abundances.tsv --report_file stats/genomes_abundances_mqc.tsv \
      --checkm_file stats/genomes_checkm_mqc.json --table_file stats/bins_general_stats_mqc.tsv

      cat ${table_header_mqc} > stats/tmp.txt && cat stats/bins_general_stats_mqc.tsv >> stats/tmp.txt \
      && mv stats/tmp.txt stats/bins_general_stats_mqc.tsv

      cat ${heatmap_header_mqc} > stats/tmp.txt && cat stats/genomes_abundances_mqc.tsv >> stats/tmp.txt \
      && mv stats/tmp.txt stats/genomes_abundances_mqc.tsv
      """
}

process ADD_QUAST_INFO_TO_BINS {
  publishDir "${params.outdir}/08_binning/08_4_mapping_on_final_bins/stats", mode: 'copy'
  
  input:
    tuple val(meta), path(bins_stat), path(quast_report)

  output:
    path "${meta.id}_bins_stat_and_quality.tsv", emit: bins_stats

  script:
    """
    awk -F"\t" '{ 
      if (NR==1) {
        val=-1; 
        for(i=1;i<=NF;i++) { 
          if (\$i ~ /${meta.id}.*/) {
            val=i;}}} 
      if(val != -1) print \$1 "\t" \$val} ' report.tsv > report_${meta.id}.tsv

    add_info_to_bin_stat.py -s $bins_stat \
                                 -q report_${meta.id}.tsv \
                                 -o "${meta.id}_bins_stat_and_quality.tsv"
    """
}

process BINS_STATS_TO_MUTLIQC_FORMAT {
  publishDir "${params.outdir}/08_binning/08_4_mapping_on_final_bins/stats", mode: 'copy'
  
  input:
    path bins_stats
    path multiqc_header_count
    path multiqc_header_size

  output:
    tuple path("bin_count_per_quality_mqc.tsv"), path("bin_size_per_quality_mqc.tsv"), emit: report

  script:
   """
   format_bins_stat_to_multiqc.py --bins_stats $bins_stats -v \
                                    --out_bins_count bin_count_per_quality.tsv \
                                    --out_bins_size bin_size_per_quality.tsv

   cat $multiqc_header_count > bin_count_per_quality_mqc.tsv
   cat bin_count_per_quality.tsv >> bin_count_per_quality_mqc.tsv

   cat $multiqc_header_size > bin_size_per_quality_mqc.tsv
   cat bin_size_per_quality.tsv >> bin_size_per_quality_mqc.tsv

   """
}
