process GENERATE_DEPTH_FILES {
  tag "${meta.id}"
  publishDir "${params.outdir}/08_binning/08_1_binning_per_sample/${meta.id}/", mode: 'copy'
  label 'BINNING'

  input:
    tuple val(meta), path(bam), path(bai)

  output:
    tuple val(meta), path("${meta.id}_depth.txt")

  script:
  """
  jgi_summarize_bam_contig_depths --outputDepth ${meta.id}_depth.txt $bam 
  """
}

process METABAT2 {
  tag "${meta.id}"
  publishDir "${params.outdir}/08_binning/08_1_binning_per_sample/${meta.id}/", mode: 'copy'
  label 'BINNING'

  input:
    tuple val(meta), path(fna), path(depth)
      
  output:
    tuple val(meta), path("metabat2/bins/"), emit: bins
    path "v_metabat2.txt", emit: v_metabat2

  script:
  """
  mkdir -p metabat2/bins/
  metabat2 --inFile $fna --abdFile $depth --outFile metabat2/bins/${meta.id}_metabat2 --numThreads ${task.cpus} --seed ${params.metabat2_seed}
  echo \$(metabat2 -h 2>&1) > v_metabat2.txt
  """
}

process MAXBIN2 {
  errorStrategy { task.exitStatus == 255 ? 'ignore' : 'finish' }
  tag "${meta.id}"
  publishDir "${params.outdir}/08_binning/08_1_binning_per_sample/${meta.id}/", mode: 'copy'
  label 'BINNING'

  input:
    tuple val(meta), path(fna), path(depth)
      
  output:
    tuple val(meta), path("maxbin2/bins/"), emit: bins
    path "v_maxbin.txt", emit: v_maxbin

  script:
  """
  cut -f1,3- $depth > ${meta.id}_maxbin2_depth.txt
  mkdir maxbin2/
  run_MaxBin.pl -contig $fna  -out maxbin2/${meta.id}_maxbin2 -abund ${meta.id}_maxbin2_depth.txt -thread ${task.cpus}

  mkdir maxbin2/bins/
  cd maxbin2/
  if [ -f ${meta.id}_maxbin2.001.fasta ]
  then 
    for file in \$(ls ${meta.id}_maxbin2*.fasta)
    do 
      mv -- "\$file" "bins/\${file%.fasta}.fa"
    done
  fi

  cd ..
  rm ${meta.id}_maxbin2_depth.txt

  run_MaxBin.pl -v &> v_maxbin.txt
  """
}

process CONCOCT {
  errorStrategy { task.exitStatus == 255 ? 'ignore' : 'finish' }
  tag "${meta.id}"
  publishDir "${params.outdir}/08_binning/08_1_binning_per_sample/${meta.id}", mode: 'copy'
  label 'BINNING'
  
  input:
    tuple val(meta), path(fna), path(bam), path(bai)
      
  output:
    tuple val(meta), path("concoct/bins/"), emit: bins
    path "v_concoct.txt", emit: v_concoct

  script:
  """
  mkdir -p concoct

  cut_up_fasta.py $fna --chunk_size 10000 --overlap_size 0 --merge_last --bedfile concoct/contigs_10K.bed > concoct/contigs_10K.fa

  concoct_coverage_table.py concoct/contigs_10K.bed  $bam > concoct/coverage_table.tsv

  concoct --composition_file concoct/contigs_10K.fa --coverage_file concoct/coverage_table.tsv --basename concoct/bins --threads ${task.cpus} 

  merge_cutup_clustering.py concoct/bins_clustering_gt1000.csv > concoct/clustering_merge.csv

  mkdir -p concoct/bins

  extract_fasta_bins.py $fna  concoct/clustering_merge.csv --output_path concoct/bins
  concoct -v &> v_concoct.txt
  
  cd concoct/bins
  for filename in *fa; 
  do 
    mv \$filename "${meta.id}_concoct_\$filename";
  done

  
  """
}

process BINETTE {
  errorStrategy { task.exitStatus == 255 ? 'ignore' : 'finish' }
  tag "${meta.id}"
  publishDir "${params.outdir}/08_binning/08_1_binning_per_sample/${meta.id}", mode: 'copy'
  label 'BINNING'

  input:
    tuple val(meta), val(bins), path(contigs)
    val min_completeness
    val max_contamination
    val(checkm2_db)
      
  output:
    tuple val(meta), path("bin_refinement/final_bins/*"), emit: bins, optional: true
    tuple val(meta), path('bins_stats.tsv'), emit: checkm_stats
    path "v_binette.txt", emit: v_binette

	script:
  """
    binette -v -m ${min_completeness} -t ${task.cpus} -d ${bins.join(' ')} -c ${contigs} -o bin_refinement/ --checkm2_db $checkm2_db

    sed -i 's/size/Size/' bin_refinement/final_bins_quality_reports.tsv 
    sed -i 's/^/${meta.id}_bin_/' bin_refinement/final_bins_quality_reports.tsv 
    sed -i 's/^${meta.id}_bin_bin_id/genome/' bin_refinement/final_bins_quality_reports.tsv 
    cut -f 2,3 --complement bin_refinement/final_bins_quality_reports.tsv > bins_stats.tsv 

    cd bin_refinement/final_bins

    if [[ \$(ls ./) ]]  #check if directoy not empty
    then 
      for filename in *fa; 
      do 
        mv \$filename "${meta.id}_\$filename"
      done
    fi

    cd ../..
    
    binette --version > v_binette.txt
  """
}

process UNBINNED_CONTIGS {
  tag "${meta.id}"
  publishDir "${params.outdir}/08_binning/08_1_binning_per_sample/${meta.id}/bin_refinement/", mode: 'copy'

  input:
    tuple val(meta), path(bins), path(assembly)
      
  output:
    path "unbinned_contigs.fasta", emit: unbinned

  """
    write_unbinned_contigs.py -b ${bins} \
    -a ${assembly} \
    -o unbinned_contigs.fasta

  """
}
