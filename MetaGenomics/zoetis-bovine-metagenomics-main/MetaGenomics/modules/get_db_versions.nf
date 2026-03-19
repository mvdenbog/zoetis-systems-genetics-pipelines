process GET_DB_VERSIONS {
  publishDir "${params.outdir}/pipeline_info", mode: 'copy'
  //label 'BINNING' // For checkm2
  label 'GENERAL'

  input:
  path host_fasta
  path kaiju_nodes
  path kraken2
  path accession
  path taxdump


  output:
    path "db_versions.tsv"

  script:
  """
  if [[ "${host_fasta}" != "" ]]
  then
    echo "Host_genome ${host_fasta}" >  host_genome_db.txt
  fi

  if [[ "${kaiju_nodes}" != "" ]] 
  then
    echo "Kaiju ${kaiju_nodes}" > kaiju_db.txt
  fi
  if [[ "${kraken2}" != "" ]] 
  then
    echo "Kaiju ${kraken2}" > kraken2_db.txt
  fi

  if [[ "${taxdump}" != "" && "${accession}" != "" ]]
  then 
    echo "Taxdump ${taxdump}" > taxdump_db.txt
    echo "Accession2taxid ${accession}" > accession_db.txt
  fi 



  if [[ `ls | grep db.txt` ]]
  then
    for i in *_db.txt
    do 
      cat \$i >> all_db.txt
    done
  else 
    touch all_db.txt
  fi
  
  db_versions.py -a all_db.txt >> db_versions.tsv
  """

}
