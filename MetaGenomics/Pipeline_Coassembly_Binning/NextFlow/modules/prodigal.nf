process PRODIGAL {
  tag "${meta.id}"

  input:
    tuple val(meta), file(assembly_file)

  output:
    tuple val(meta), path("prodigal.gff"), emit: gff
    tuple val(meta), path("prodigal.faa"), emit: faa
    path "v_prodigal.txt", emit: v_prodigal

  script:
  """
  prodigal -i ${assembly_file} -c -p meta -f gff -a prodigal.faa -o prodigal.gff 
  
  prodigal -v 2> v_prodigal.txt

  """
}