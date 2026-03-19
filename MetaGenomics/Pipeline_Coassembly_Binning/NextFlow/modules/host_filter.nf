process HOST_FILTER {
  tag "${meta.id}"

  publishDir "${params.outdir}/01_clean_qc/01_1_cleaned_reads/", mode: 'copy', pattern: 'cleaned_*.fastq.gz'
  publishDir "${params.outdir}/01_clean_qc/01_1_cleaned_reads/", mode: 'copy', pattern: '*.bam'
  publishDir "${params.outdir}/01_clean_qc/01_1_cleaned_reads/logs", mode: 'copy',
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

process HOST_FILTER_HIFI {
  tag "${meta.id}"
  label "MINIMAP2"

  publishDir "${params.outdir}/01_clean_qc/01_1_cleaned_reads/", mode: 'copy', pattern: 'cleaned_*.fastq.gz'
  publishDir "${params.outdir}/01_clean_qc/01_1_cleaned_reads/", mode: 'copy', pattern: '*.bam'
  publishDir "${params.outdir}/01_clean_qc/01_1_cleaned_reads/logs", mode: 'copy',
    saveAs: {filename -> 
    if (filename.indexOf(".flagstat") > 0) "$filename"
    else null}

  input:
    tuple val(meta), path(reads)
    path fasta

  output:
    tuple val(meta), path("cleaned_${meta.id}.fastq.gz"), emit: reads
    path "host_filter_flagstat/${meta.id}.host_filter.flagstat", emit: hf_report
    path "${meta.id}.no_filter.flagstat", emit: nf_report
    path "v_minimap2.txt", emit: v_minimap
    path "v_samtools.txt", emit: v_samtools


  script:
    """
    minimap2 -ax map-hifi -t ${task.cpus} ${fasta} ${reads} | samtools sort  -@ ${task.cpus} -o ${meta.id}.bam
    
    samtools view -@ ${task.cpus} -bh -f 4 ${meta.id}.bam > ${meta.id}.without_host.bam
    
    mkdir host_filter_flagstat
    
    samtools flagstat ${meta.id}.bam -@ ${task.cpus} > ${meta.id}.no_filter.flagstat
    samtools flagstat ${meta.id}.without_host.bam -@ ${task.cpus} > host_filter_flagstat/${meta.id}.host_filter.flagstat
    samtools fastq  -@ ${task.cpus} ${meta.id}.without_host.bam | gzip > cleaned_${meta.id}.fastq.gz
    
    rm ${meta.id}.bam
    rm ${meta.id}.without_host.bam

    minimap2 --version &> v_minimap2.txt
    samtools --version &> v_samtools.txt
    """
}
