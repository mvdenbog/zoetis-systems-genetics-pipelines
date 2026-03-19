process CHUNK_ASSEMBLY_FILTER {
  label 'ASSEMBLY_FILTER'

  publishDir "${params.outdir}/03_filtering/chunk_assembly_filter/", mode: 'copy'  
  input:
    tuple val(meta), path(assembly_file), path(idxstats)
    val min_cpm

  output:
    tuple val(meta), path("${chunk_name}_select_cpm${min_cpm}.fasta"), emit: chunk_selected
    tuple val(meta), path("${chunk_name}_discard_cpm${min_cpm}.fasta"), emit: chunk_discarded

  script:
    chunk_name = "${meta.id}_${assembly_file.baseName}" // assembly_file.baseName
    """
    filter_contig_per_cpm.py -v -i ${idxstats} -f ${assembly_file} -c ${min_cpm} -s ${chunk_name}_select_cpm${min_cpm}.fasta -d ${chunk_name}_discard_cpm${min_cpm}.fasta
    """
}

process MERGE_ASSEMBLY_FILTER {
  label 'ASSEMBLY_FILTER'

  tag "${meta.id}"
  publishDir "${params.outdir}/${publishDir_path}/", mode: 'copy', pattern: "*_select_contigs*"
  publishDir "${params.outdir}/${publishDir_path}/discard_contigs", mode: 'copy', pattern: "*_discard_contigs*"

  input:
    tuple val(meta), path(select_fasta)
    tuple val(meta), path(discard_fasta)
    val min_cpm
    val(publishDir_path)

  output:
    tuple val(meta), path("${meta.id}_select_contigs_cpm${min_cpm}.fasta"), emit: merged_selected
    tuple val(meta), path("${meta.id}_discard_contigs_cpm${min_cpm}.fasta"), emit: merged_discarded

  shell:
    '''
    echo !{select_fasta} | sed "s/ /\\n/g" | sort > select_list
    echo !{discard_fasta} | sed "s/ /\\n/g" | sort > discard_list

    for i in `cat select_list` ; 
    do 
      cat $i >> !{meta.id}_select_contigs_cpm!{min_cpm}.fasta
    done
    
    for j in `cat discard_list` ; 
    do 
      cat $j >> !{meta.id}_discard_contigs_cpm!{min_cpm}.fasta
    done
    
    rm select_list
    rm discard_list
    '''
}
