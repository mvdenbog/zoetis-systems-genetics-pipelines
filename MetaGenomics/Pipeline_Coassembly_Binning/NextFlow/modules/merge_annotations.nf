process MERGE_ANNOTATIONS {
  publishDir "${params.outdir}/04_structural_annot/${meta.id}/", mode: 'copy'
  tag "${meta.id}"

  input:
    tuple val(meta), file(assembly_file), file(faa_file), file(cds_gff), file(rrna_gff), file(trna_gff)

  output:
    tuple val(meta), file("${meta.id}.gff"), emit: gff
    tuple val(meta), file("${meta.id}.ffn"), emit: ffn
    tuple val(meta), file("${meta.id}.faa"), emit: faa
    path "${meta.id}.txt", emit: report

  script:
  """
    merge_annotations.py -c $cds_gff -r $rrna_gff -t $trna_gff -v  \
                         --contig_seq $assembly_file --faa_file $faa_file \
                         --ffn_output ${meta.id}.ffn  --gff_output ${meta.id}.gff --faa_output ${meta.id}.faa \
                         --report_output ${meta.id}.txt
  """
}