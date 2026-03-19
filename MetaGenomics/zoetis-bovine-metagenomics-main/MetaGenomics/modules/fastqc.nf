process FASTQC {
  tag "${meta.id}"
  //label 'FASTQC'
  label 'GENERAL'
//container 'europe-west1-docker.pkg.dev/zoetis-bovine-metagenomics/containers/metagwgs:1.5'

  publishDir "${params.outdir}/${meta.id}/CleanQC/", mode: 'copy', pattern : "$outdir/${meta.id}/*"
  input:
    tuple val(meta), path(reads)

  output:
    path "$outdir/${meta.id}/*.zip", emit: zip
    path "$outdir/${meta.id}/*.html", emit: html
    path "v_fastqc.txt", emit: v_fastqc 

  script:
    if (reads[0]==~/cleaned.*/){ outdir="fastqc_cleaned" } else { outdir="fastqc_raw" }
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
}
