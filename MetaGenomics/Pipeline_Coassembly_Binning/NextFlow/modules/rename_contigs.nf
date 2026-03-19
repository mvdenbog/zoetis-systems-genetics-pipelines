process RENAME_CONTIGS {
  tag "${meta.id}"
  publishDir "${params.outdir}/02_assembly/02_1_primary_assembly", mode: 'copy'
  label 'PYTHON'

  input:
    tuple val(meta), path("${meta.id}.raw.fna")
    


  output:
    tuple val(meta), path("${meta.id}.fna"), emit: fna
    path("${meta.id}_original_to_new_contig_name.tsv")

  script:
  """
    rename_contigs.py --sample ${meta.id}  --fna_file ${meta.id}.raw.fna --out_fna ${meta.id}.fna -v -t ${meta.id}_original_to_new_contig_name.tsv 

  """
}