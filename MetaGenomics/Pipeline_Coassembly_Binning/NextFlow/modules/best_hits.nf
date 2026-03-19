process BEST_HITS {
  publishDir "${params.outdir}/06_func_annot/06_3_functional_annotation", mode: 'copy'

  input:
    tuple val(meta), path(m8)

  output:
    path "${meta.id}.best_hit", emit: best_hits

  script:
    """
    filter_diamond_hits.py -o ${meta.id}.best_hit ${m8}
    """
}