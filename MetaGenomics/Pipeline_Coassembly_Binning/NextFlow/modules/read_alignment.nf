process BWA_MEM {
   tag "${meta.id}"
   publishDir "${params.outdir}/$publishDir_path/${meta.id}", mode: 'copy', pattern:"*bam*"

   input:
      tuple val(meta), path(fna), path(reads)
      val(publishDir_path)


   output:
      tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
      path "v_bwa.txt", emit: v_bwa_mem2
      path "v_samtools.txt", emit: v_samtools

   script:
      """
      bwa-mem2 index ${fna} -p ${fna}
      bwa-mem2 mem -t ${task.cpus} ${fna} ${reads[0]} ${reads[1]} | samtools view -@ ${task.cpus} -bS - | samtools sort -@ ${task.cpus} - -o ${meta.id}.bam
      
      samtools index -@ ${task.cpus} ${meta.id}.bam

      bwa-mem2 version > v_bwa.txt
      samtools --version &> v_samtools.txt
      """
}

process MINIMAP2 {
   tag "${meta.id}"
   label 'MINIMAP2'
   publishDir "${params.outdir}/$publishDir_path/${meta.id}", mode: 'copy', pattern:"*bam*"

   input:
      tuple val(meta), path(fna), path(reads)
      val(publishDir_path)

   output:
      tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
      path "v_minimap2.txt", emit: v_minimap2
      path "v_samtools.txt", emit: v_samtools

   script:
      """
      # align reads to contigs, keep only primary aln and sort resulting bam 
      minimap2 -t ${task.cpus} -ax map-hifi $fna $reads | samtools view -@ ${task.cpus} -b -F 2304 | samtools sort -@ ${task.cpus} -o ${meta.id}.bam  
      
      samtools index ${meta.id}.bam -@ ${task.cpus}

      samtools --version &> v_samtools.txt
      minimap2 --version &> v_minimap2.txt
      """
}

