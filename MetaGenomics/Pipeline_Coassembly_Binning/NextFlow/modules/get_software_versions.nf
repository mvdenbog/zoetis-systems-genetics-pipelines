process GET_SOFTWARE_VERSIONS {
  publishDir "${params.outdir}/pipeline_info", mode: 'copy',
  saveAs: {filename ->
    if (filename.indexOf(".csv") > 0) filename
    else null
  }

  input:
    path software_versions

  output:
    path 'software_versions_mqc.yaml', emit: yaml
    path "software_versions.csv"

  script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    python --version &> v_python.txt
    multiqc --version &> v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}