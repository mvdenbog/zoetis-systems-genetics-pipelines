process FASTQC {
  tag "${meta.id}"
  label 'FASTQC'

  publishDir "${params.outdir}/01_clean_qc/01_2_qc/", mode: 'copy', pattern : "$outdir/${meta.id}/*"
  input:
    tuple val(meta), path(reads)

  output:
    path "$outdir/${meta.id}/*.zip", emit: zip
    path "$outdir/${meta.id}/*.html", emit: html
    path "v_fastqc.txt", emit: v_fastqc 

  script:
    if (reads[0]==~/cleaned.*/){ outdir="fastqc_cleaned" } else { outdir="fastqc_raw" }
    if (meta.type=="SR"){ 
      option = "--nogroup" 
      """
      if [[ "$outdir" == "fastqc_raw" && ! -f ${meta.id}_R1.fastq.gz ]]
      then
        ln -s ${reads[0]} ${meta.id}_R1.fastq.gz
        ln -s ${reads[0]} ${meta.id}_R2.fastq.gz
        fastq="${meta.id}_R1.fastq.gz ${meta.id}_R2.fastq.gz"
      else 
        fastq="${reads[0]} ${reads[1]}"
      fi

      mkdir -p $outdir/${meta.id} 
      fastqc --nogroup --quiet -o $outdir/${meta.id} --threads ${task.cpus} \$fastq

      fastqc --version &> v_fastqc.txt
      """
    } else if (meta.type=="HIFI"){ 
      """
      if [[ "$outdir" == "fastqc_raw" && ! -f ${meta.id}.fastq.gz ]]
      then
        ln -s ${reads} ${meta.id}.fastq.gz
        fastq="${meta.id}.fastq.gz"
      else 
        fastq="${reads}"
      fi

      mkdir -p $outdir/${meta.id} 
      fastqc --quiet -o $outdir/${meta.id} --threads ${task.cpus} \$fastq

      fastqc --version &> v_fastqc.txt
      """
    }
}
