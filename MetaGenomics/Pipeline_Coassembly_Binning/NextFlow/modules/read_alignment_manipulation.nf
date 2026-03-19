process GET_ALIGNMENT_METRICS {
   tag "${meta.id}"
   publishDir "${params.outdir}/$publishDir_path/${meta.id}", mode: 'copy'

   input:
      tuple val(meta), path(bam), path(bai)
      val(publishDir_path)

   output:
      tuple val(meta), path("${meta.id}.coverage.tsv"), emit: sam_coverage
      tuple val(meta), path("${meta.id}.idxstats"), emit: sam_idxstat
      path "${meta.id}.flagstat", emit: sam_flagstat
      path "${meta.id}*"
      path "v_samtools.txt", emit: v_samtools

   script:
      """
      samtools flagstat -@ ${task.cpus}  ${bam} > ${meta.id}.flagstat
      samtools coverage  ${bam} > ${meta.id}.coverage.tsv
      samtools idxstats ${bam} > ${meta.id}.idxstats

      samtools --version &> v_samtools.txt
      """
}

