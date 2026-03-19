process BARRNAP {
  tag "${meta.id}"

  input:
    tuple val(meta), file(assembly_file)

  output:
    tuple val(meta), path("barrnap.gff"), emit: gff
    path "v_barrnap.txt", emit: v_barrnap

  script:
  """
    barrnap --threads ${task.cpus} ${assembly_file} > barrnap.gff

    barrnap --version  2> v_barrnap.txt

  """
}
