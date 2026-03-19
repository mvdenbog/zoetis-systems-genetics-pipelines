process MERGE_QUANT_ANNOT_BEST {
  publishDir "${params.outdir}/06_func_annot/06_3_functional_annotation", mode: 'copy'

  input:
    path quant_table
    path annot
    path best_hits

  output:
    path "Quantifications_and_functional_annotations.tsv", emit: merged

  script:
    """
    awk '{
        if(NR == 1) {
            print \$0 "\t" "sum"}
        else {
            for (i=1; i<=NF; i++)
            {
                if (i == 1)
                {
                    sum = O;
                }
                else {
                    sum = sum + \$i;
                }
            }
            print \$0 "\t" sum
        }
    }' ${quant_table} > ${quant_table}.sum
    ls ${annot} | cat > List_of_functionnal_annotations_files.txt
    ls ${best_hits} | cat > List_of_diamond_files.txt
    merge_abundance_and_functional_annotations.py -t ${quant_table}.sum -f List_of_functionnal_annotations_files.txt -d List_of_diamond_files.txt -o Quantifications_and_functional_annotations.tsv
    """
}