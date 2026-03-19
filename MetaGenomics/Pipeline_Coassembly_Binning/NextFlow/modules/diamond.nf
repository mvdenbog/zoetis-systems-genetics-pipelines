process DIAMOND {
    publishDir "${params.outdir}/05_protein_alignment/05_1_database_alignment/$meta.id", mode: 'copy', pattern: "*.m8"
    tag "${meta.id}"

    input:
      tuple val(meta), path(faa)
      path diamond_bank

    output:
      tuple val(meta), path("${meta.id}_aln_diamond.m8"), emit: m8
      path "v_diamond.txt", emit: v_diamond

    script:
      fmt="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle"
      fmt_tab=fmt.replaceAll(" ","\t")
      """
      echo "$fmt_tab" > head.m8
      diamond blastp -p ${task.cpus} -d ${diamond_bank} -q ${faa} -o ${meta.id}_aln_diamond.nohead.m8 -f 6 $fmt
      cat head.m8 ${meta.id}_aln_diamond.nohead.m8 > ${meta.id}_aln_diamond.m8
      rm ${meta.id}_aln_diamond.nohead.m8
      rm head.m8

      diamond help &> v_diamond.txt
      """
}
