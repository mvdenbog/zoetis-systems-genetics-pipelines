process ASSIGN_TAXONOMY {
  tag "${meta.id}"
  publishDir "${params.outdir}/07_taxo_affi/07_1_affiliation_per_sample/${meta.id}", mode: 'copy'
  label 'PYTHON'

  input:
    tuple path(accession2taxid), path(new_taxdump)
    tuple val(meta), path(m8), path(sam_coverage)

  output:
    tuple val(meta.id), path("${meta.id}.percontig.tsv"), emit: t_percontig
    tuple val(meta.id), path("${meta.id}.pergene.tsv"), emit: t_pergene
    tuple val(meta.id), path("${meta.id}.warn.tsv"), emit: t_warn
    path "${meta.id}_quantif_percontig.tsv", emit: q_all
    path "${meta.id}_quantif_percontig_by_superkingdom.tsv", emit: q_superkingdom
    path "${meta.id}_quantif_percontig_by_phylum.tsv", emit: q_phylum
    path "${meta.id}_quantif_percontig_by_order.tsv", emit: q_order
    path "${meta.id}_quantif_percontig_by_class.tsv", emit: q_class
    path "${meta.id}_quantif_percontig_by_family.tsv", emit: q_family
    path "${meta.id}_quantif_percontig_by_genus.tsv", emit: q_genus
    path "${meta.id}_quantif_percontig_by_species.tsv", emit: q_species
    path "top_taxons_per_contig.tsv", emit: top_taxon_file
    path "krona_mean_depth_abundance/${meta.id}.krona", emit: krona_depth
    path "krona_reads_count_abundance/${meta.id}.krona", emit: krona_reads_count

  script:
    """
    new_taxdump_var=$new_taxdump
    
    if [ "\${new_taxdump_var#*.}" == "tar.gz" ]
    then 
    
     	echo "$new_taxdump taxdump dir is tar.gz archive."
    	mkdir new_taxdump
     	tar xzvf $new_taxdump merged.dmp nodes.dmp taxidlineage.dmp names.dmp
     	mv merged.dmp nodes.dmp taxidlineage.dmp names.dmp new_taxdump/
     	new_taxdump_var="new_taxdump"
    fi


    aln_to_tax_affi.py -a ${accession2taxid} --taxonomy \$new_taxdump_var \
                  -o ${meta.id} -b ${m8} --keep_only_best_aln  \
                   -v --write_top_taxons 
    merge_contig_quantif_perlineage.py -c ${meta.id}.percontig.tsv -s ${sam_coverage} -o ${meta.id}
    
    new_taxdump_original=$new_taxdump
    if [ "\${new_taxdump_original#*.}" == "tar.gz" ]
    then 
    	rm -r new_taxdump
    fi 
    
    """
}


process  PLOT_TAXONOMIC_AFFILIATIONS{
  publishDir "${params.outdir}/07_taxo_affi/07_3_plot/", mode: 'copy'
  label 'PYTHON'

  input:
    // path "*_quantif_percontig.tsv"
    path quantif_percontig

  output:
    path "*html"

  script:
    """
    plot_contigs_taxonomic_affiliation.py $quantif_percontig --output_dir "."


    """
}


