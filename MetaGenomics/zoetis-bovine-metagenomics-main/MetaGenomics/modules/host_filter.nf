process HOST_FILTER {
  tag "${meta.id}"
  label 'GENERAL'
  publishDir "${params.outdir}/${meta.id}/CleanedReads/", mode: 'copy', pattern: 'cleaned_*.fastq.gz'
  publishDir "${params.outdir}/${meta.id}/CleanedReads/", mode: 'copy', pattern: '*.bam'
  publishDir "${params.outdir}/${meta.id}/CleanedReads/logs", mode: 'copy',
    saveAs: {filename -> 
    if (filename.indexOf(".flagstat") > 0) "$filename"
    else null}

  input:
    tuple val(meta), path(reads)
    path fasta
    path index

  output:
    tuple val(meta), path("cleaned_${meta.id}*.fastq.gz"), emit: reads
    path "host_filter_flagstat/${meta.id}.host_filter.flagstat", emit: hf_report
    path "${meta.id}.no_filter.flagstat", emit: nf_report
    path "v_bwa.txt", emit: v_bwa_mem2
    path "v_samtools.txt", emit: v_samtools

  script:
    """
    bwa-mem2 mem -t ${task.cpus} ${fasta} ${reads[0]} ${reads[1]} > ${meta.id}.bam
    samtools view -bhS -f 12 ${meta.id}.bam > ${meta.id}.without_host.bam
    mkdir host_filter_flagstat
    samtools flagstat ${meta.id}.bam > ${meta.id}.no_filter.flagstat
    samtools flagstat ${meta.id}.without_host.bam >> host_filter_flagstat/${meta.id}.host_filter.flagstat
    samtools sort -n -o ${meta.id}.without_host_sort.bam ${meta.id}.without_host.bam     
    samtools fastq -N -1 cleaned_${meta.id}_R1.fastq.gz -2 cleaned_${meta.id}_R2.fastq.gz ${meta.id}.without_host_sort.bam
    rm ${meta.id}.bam
    rm ${meta.id}.without_host.bam
    rm ${meta.id}.without_host_sort.bam


    bwa-mem2 version > v_bwa.txt
    samtools --version &> v_samtools.txt
    """
}

