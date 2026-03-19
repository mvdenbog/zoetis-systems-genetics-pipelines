process EGGNOG_MAPPER {
  publishDir "${params.outdir}/06_func_annot/06_3_functional_annotation", mode: 'copy', pattern: "${meta.id}_diamond_one2one.emapper.*"
  tag "${meta.id}"
  label 'EGGNOG'

  input:
  tuple val(meta), path(faa)
  path db

  output:
  path "${meta.id}_diamond_one2one.emapper.seed_orthologs", emit: seed
  path "${meta.id}_diamond_one2one.emapper.annotations", emit: annot
  path 'v_eggnogmapper.txt', emit: version

  script:
  """
  emapper.py -i ${faa} --output ${meta.id}_diamond_one2one -m diamond --cpu ${task.cpus} --data_dir ${db} --target_orthologs one2one
  emapper.py -v &> v_eggnogmapper.txt
  """
}