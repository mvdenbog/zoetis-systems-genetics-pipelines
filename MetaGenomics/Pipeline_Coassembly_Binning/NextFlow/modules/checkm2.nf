process CHECKM2 {
  label 'BINNING'
  publishDir "${params.outdir}/08_binning/08_1_binning_per_sample/${meta.id}", mode: 'copy', pattern: "checkm2/*"
  input:
    tuple val(meta), val(bins)
    val(checkm2_db)
      
  output:
    tuple val(meta.id), path("checkm2/*"), emit: checkm2
    path "v_checkm2.txt", emit: v_checkm2

  script:
  """
    checkm2 predict -x fa --threads ${task.cpus} --database_path ${checkm2_db} --input ${bins.join(' ')} --output-directory checkm2/ 

    echo "\$(checkm2 --version)_modified" > v_checkm2.txt
  """ 
}