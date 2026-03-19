process CIRCULAR_CONTIGS_METAFLYE {
  tag "${meta.id}"
  publishDir "${params.outdir}/02_assembly/02_1_primary_assembly", mode: 'copy'
  label 'PYTHON'

  input:
    tuple val(meta), path("${meta.id}.fna"), path(infos)

  output:
    tuple val(meta), path("circular_contigs/"), emit: circular

  script:
    """
    retrieve_circular_contigs.py -a 'metaflye' -f ${meta.id}.fna -i ${infos} -o circular_contigs/
    """
}


process CIRCULAR_CONTIGS_HIFIASM {
  tag "${meta.id}"
  publishDir "${params.outdir}/02_assembly/02_1_primary_assembly", mode: 'copy'
  label 'PYTHON'

  input:
    tuple val(meta), path("${meta.id}.fna")

  output:
    tuple val(meta), path("circular_contigs/"), emit: circular

  script:
    """
    retrieve_circular_contigs.py -a 'hifiasm-meta' -f ${meta.id}.fna -o circular_contigs/
    """
}