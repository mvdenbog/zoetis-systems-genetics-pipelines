taxo_list = "all superkingdom phylum class order family genus species"

process QUANTIF_AND_TAXONOMIC_TABLE_CONTIGS {
  publishDir "${params.outdir}/07_taxo_affi/07_2_affiliation_merged", mode: 'copy'
  label 'PYTHON'

  input:
    path q_all
    path q_superkingdom
    path q_phylum
    path q_order
    path q_class
    path q_family
    path q_genus
    path q_species

  output:
    path "quantification_by_contig_lineage*.tsv", emit: quantif_by_contig_lineage

  script:
    """
      echo "${q_all}" > all.txt
      echo "${q_superkingdom}" > superkingdom.txt
      echo "${q_phylum}" > phylum.txt
      echo "${q_order}" > order.txt
      echo "${q_class}" > class.txt
      echo "${q_family}" > family.txt
      echo "${q_genus}" > genus.txt
      echo "${q_species}" > species.txt
      for i in ${taxo_list} ;
      do
      quantification_by_contig_lineage.py -i \$i".txt" -o quantification_by_contig_lineage_\$i".tsv"
      done
    """
}