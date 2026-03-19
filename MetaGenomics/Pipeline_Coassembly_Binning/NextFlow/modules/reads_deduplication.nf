process READS_DEDUPLICATION {
  tag "${meta.id}"
  publishDir "${params.outdir}/02_assembly/02_2_deduplicated_reads/", mode: 'copy', pattern: '*.fastq.gz'
  publishDir "${params.outdir}/02_assembly/02_3_reads_vs_primary_assembly/${meta.id}/", mode: 'copy', pattern: '*.bam*'
  
  input:
  tuple val(meta), path(assembly), path(reads)

  output:
  tuple val(meta), path("${meta.id}*_dedup.fastq.gz"), emit: dedup_reads
  tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
  
  path "${meta.id}_singletons.fastq.gz", emit: singletons

  path "v_bwa.txt", emit: v_bwa_mem2
  path "v_samtools.txt", emit: v_samtools

  script:
  """
  # Align reads against assembly
  bwa-mem2 index ${assembly} -p ${assembly}
  bwa-mem2 mem -t ${task.cpus}  ${assembly} ${reads[0]} ${reads[1]} | samtools view -@ ${task.cpus} -bS - | samtools sort -@ ${task.cpus} -n -o ${meta.id}.sort.bam -
  
  # Identify and removed duplicated reads from the bam file
  samtools fixmate -@ ${task.cpus} -m ${meta.id}.sort.bam ${meta.id}.fixmate.bam
  samtools sort -@ ${task.cpus} -o ${meta.id}.fixmate.positionsort.bam ${meta.id}.fixmate.bam
  samtools markdup -@ ${task.cpus}  -r -S -s -f ${meta.id}.stats ${meta.id}.fixmate.positionsort.bam ${meta.id}.filtered.bam
  
  # final bam file without duplicated reads
  samtools sort -@ ${task.cpus} -o ${meta.id}.bam ${meta.id}.filtered.bam
  samtools index -@ ${task.cpus}  ${meta.id}.bam 

  # Get deduplicated reads
  samtools sort -@ ${task.cpus}  -n -o ${meta.id}.filtered.n_sorted.bam ${meta.id}.bam  
  samtools fastq -@ ${task.cpus}  -N -s ${meta.id}_singletons.fastq.gz -1 ${meta.id}_R1_dedup.fastq.gz -2 ${meta.id}_R2_dedup.fastq.gz ${meta.id}.filtered.n_sorted.bam 
  
  # clean directory
  rm ${meta.id}.sort.bam
  rm ${meta.id}.fixmate.bam
  rm ${meta.id}.fixmate.positionsort.bam
  rm ${meta.id}.filtered.n_sorted.bam
  rm ${meta.id}.filtered.bam


  bwa-mem2 version > v_bwa.txt
  samtools --version &> v_samtools.txt
  """
}