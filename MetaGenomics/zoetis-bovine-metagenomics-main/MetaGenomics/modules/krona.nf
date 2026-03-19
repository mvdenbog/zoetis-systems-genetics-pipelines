process PLOT_KRONA {

  publishDir "${params.outdir}/$publishDir_path", mode: 'copy', pattern: "*html"
  input:
    path krona_tab_files
    val(publishDir_path)
    val(plot_name)

  output:
    path plot_name
    path "v_kronatools.txt", emit: v_kronatools


  script:
    """
    ktImportText -o $plot_name $krona_tab_files

    ktImportText &> v_kronatools.txt
    """
}