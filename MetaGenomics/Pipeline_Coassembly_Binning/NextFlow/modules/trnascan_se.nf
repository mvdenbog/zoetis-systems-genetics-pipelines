process TRNASCAN_SE {
  tag "${meta.id}"

  input:
    tuple val(meta), file(assembly_file)

  output:
    tuple val(meta), path("trnascan_se.gff"), emit: gff
    path "v_tRNAscan.txt", emit: v_tRNAscan

  script:
  """
    tRNAscan-SE -Q -B --gff trnascan_se.gff  --thread ${task.cpus} --stats trnascan_se.log  ${assembly_file}

    tRNAscan-SE -h 2> v_tRNAscan.txt


  """
}
