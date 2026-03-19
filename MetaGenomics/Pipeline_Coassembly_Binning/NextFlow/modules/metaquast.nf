process QUAST {
    label 'QUAST'
    publishDir "${params.outdir}/${outdir}/", mode: 'copy'

    input:
        path(assemblies)
        val outdir

    output:
        path("assembly_metric/*"), emit: all
        path("assembly_metric/report.tsv"), emit: report
        path("v_quast.txt"), emit: v_quast

    script:
        """
        metaquast.py --threads ${task.cpus} --rna-finding --max-ref-number 0 --min-contig 0 ${assemblies} -o assembly_metric/

        quast -v &> v_quast.txt
        """
}
