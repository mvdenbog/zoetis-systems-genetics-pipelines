process MULTIQC {
  publishDir "${params.outdir}/MultiQC", mode: 'copy'

  input:
    path multiqc_config
    path 'software_versions/*'
    path 'cutadapt_report/'
    path 'sickle_report/'
    path 'reads_on_host_genome/'
    path 'fastqc_raw_report/'
    path 'fastqc_clean_report/'
    path 'kaiju_report/'
    path 'kraken2_report/'
    path "unfiltered_assembly_flagstat/"

  output:
    path "multiqc_report.html", emit: report
    path "multiqc_data/*"
  
  script:
    """
    multiqc . --config ${multiqc_config} -m custom_content -m fastqc -m cutadapt -m sickle -m kaiju -m quast -m prokka -m featureCounts -m samtools -m kraken
    """
}

