// Quantification of reads on each gene in each sample.
process FEATURE_COUNTS {
  tag "${meta.id}"
  label 'QUANTIFICATION'
  publishDir "${params.outdir}/06_func_annot/06_2_quantification", mode: 'copy', pattern: "${meta.id}.featureCounts*"

  input:
    tuple val(meta), file(gff_prokka), file(bam), file(bam_index)

  output:
    path "${meta.id}.featureCounts.tsv", emit: count_table
    path "${meta.id}.featureCounts.tsv.summary", emit: summary
    path "${meta.id}.featureCounts.stdout"
    path "v_featurecounts.txt", emit: v_featurecounts

  script:
    if (meta.type=="SR"){ option = "-p --countReadPairs" } else { option = "-L" }
    """
    featureCounts -T ${task.cpus} $option -O -t gene -g ID -a ${gff_prokka} -o ${meta.id}.featureCounts.tsv ${bam} &> ${meta.id}.featureCounts.stdout
    featureCounts -v &> v_featurecounts.txt
    """
}

// Create table with sum of reads for each global cluster of genes in each sample.
process QUANTIFICATION_TABLE {
  label 'PYTHON'

  input:
    path clusters_contigs
    path global_clusters_clusters
    path counts_files

  output:
    path "Clusters_Count_table_all_samples.txt", emit: quantification_table
    path "Correspondence_global_clstr_genes.txt"

  script:
    """
    ls ${clusters_contigs} | cat > List_of_contigs_files.txt
    ls ${counts_files} | cat > List_of_count_files.txt
    quantification_clusters.py -t ${global_clusters_clusters} -l List_of_contigs_files.txt -c List_of_count_files.txt -oc Clusters_Count_table_all_samples.txt -oid Correspondence_global_clstr_genes.txt
    """
}

workflow QUANTIFICATION {
 
  take:
    ch_gff_and_bam  // channel: [ val(meta), path(gff), path(bam), path(bam_index) ]
    ch_individual_clstr_table
    ch_global_clstr_table

  main:
  
    FEATURE_COUNTS(ch_gff_and_bam)
    ch_count_table = FEATURE_COUNTS.out.count_table.collect()
    ch_quant_report = FEATURE_COUNTS.out.summary
    QUANTIFICATION_TABLE(ch_individual_clstr_table.collect(), ch_global_clstr_table.collect(), ch_count_table)

  emit:
    quantification_table = QUANTIFICATION_TABLE.out.quantification_table
    quant_report = ch_quant_report
    v_featurecounts = FEATURE_COUNTS.out.v_featurecounts
}


