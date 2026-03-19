process FUNCTIONAL_ANNOT_TABLE {
  publishDir "${params.outdir}/06_func_annot/06_3_functional_annotation", mode: 'copy'

  input:
    path merged_quant_annot_best

  output:
    path "*", emit: functional_annot

  script:
    """
    quantification_by_functional_annotation.py -i ${merged_quant_annot_best}
    """
}