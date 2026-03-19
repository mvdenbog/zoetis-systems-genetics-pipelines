process SICKLE {
  tag "${meta.id}"

  publishDir "${params.outdir}/01_clean_qc/01_1_cleaned_reads/", mode: 'copy', pattern: 'cleaned_*.fastq.gz'
  publishDir "${params.outdir}/01_clean_qc/01_1_cleaned_reads/logs", mode: 'copy', pattern: '*_sickle.log'

  input:
    tuple val(meta), path(reads)

  output:
    tuple val(meta), path("*${meta.id}*_R*.fastq.gz"), emit: reads
    path "${meta.id}_single_sickle.fastq.gz", emit: single
    path "${meta.id}_sickle.log", emit: report
    path "v_sickle.txt", emit : v_sickle

  script:


    if (params.skip_host_filter) {
      // output are final cleaned files
      options = "-o cleaned_${meta.id}_R1.fastq.gz -p cleaned_${meta.id}_R2.fastq.gz"
    }
    else {
      //tempory files not saved in publish dir
      options = "-o ${meta.id}_sickle_R1.fastq.gz -p ${meta.id}_sickle_R2.fastq.gz"
    }
    options += " -t " + params.quality_type
    """
    sickle 'pe' -f ${reads[0]} -r ${reads[1]} $options \
    -s ${meta.id}_single_sickle.fastq.gz -g > ${meta.id}_sickle.log

    sickle --version &> v_sickle.txt
    """
}