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
    path "unfiltered_assembly_flagstat/"
    path 'quast_primary/report.tsv'
    path 'quast_filtered/report.tsv'
    path "final_assembly_flagstat/"
    path 'prokka_report/'
    path "featureCounts_report/"
    path 'binning_stat/'
    path 'binning_stat/'

  output:
    path "multiqc_report.html", emit: report
    path "multiqc_data/*"
  
  script:
    """
    multiqc . --config ${multiqc_config} -m custom_content -m fastqc -m cutadapt -m sickle -m kaiju -m quast -m prokka -m featureCounts -m samtools
    """
}


