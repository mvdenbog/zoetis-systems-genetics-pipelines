process GET_DB_VERSIONS {
  publishDir "${params.outdir}/pipeline_info", mode: 'copy'
  label 'BINNING' // For checkm2

  input:
  path host_fasta
  path kaiju_nodes
  path db_eggnog_mapper
  path accession
  path taxdump
  path diamond
  path gtdbtk
  path checkm2


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

  if [[ "${db_eggnog_mapper}" != "" ]]
  then   
    echo "Eggnog_Mapper ${db_eggnog_mapper}" > eggnog_db.txt
  fi

  if [[ "${taxdump}" != "" && "${accession}" != "" ]]
  then 
    echo "Taxdump ${taxdump}" > taxdump_db.txt
    echo "Accession2taxid ${accession}" > accession_db.txt
  fi 

  if [[ "${diamond}" != "" ]]
  then   
    echo "Diamond ${diamond}" > diamond_db.txt
  fi

  if [[ "${gtdbtk}" != "" ]]
  then   
    echo "GTDBTK ${gtdbtk}" > gtdbtk_db.txt
  fi

  if [[ "${checkm2}" != "" ]]
  then   
    echo "Checkm2 ${checkm2}" > checkm2_db.txt
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
