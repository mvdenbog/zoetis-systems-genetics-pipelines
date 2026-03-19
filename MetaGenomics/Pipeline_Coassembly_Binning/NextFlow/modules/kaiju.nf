include {PLOT_KRONA} from  './krona'

taxon_levels = "phylum class order family genus species"

process KAIJU {
  tag "${meta.id}"
  label "KAIJU"

  publishDir "${params.outdir}/01_clean_qc/01_3_taxonomic_affiliation_reads", mode: 'copy', pattern: '*_kaiju.out.gz'


  input:
    tuple val(meta), path(reads)
    tuple path(nodes), path(fmi), path(names)

  output:
    path  "${meta.id}_kaiju.out.gz",  emit: kaiju_result
    path "${meta.id}.krona", emit: krona_tab_file
    path "*.summary_*", emit: k_all
    path "*.summary_species", emit: k_species
    path "*.summary_genus", emit: k_genus
    path "*.summary_family", emit: k_family
    path "*.summary_class", emit: k_class
    path "*.summary_order", emit: k_order
    path "*.summary_phylum", emit: k_phylum
    path "v_kaiju.txt", emit: v_kaiju

  script:
    if (meta.type=="SR"){ 
      input = " -i ${reads[0]} -j ${reads[1]}"
    } else if (meta.type=="HIFI"){ 
      input = "-i ${reads}"
    }
    """
    kaiju -z ${task.cpus} -t ${nodes} -f ${fmi} $input -o ${meta.id}_kaiju_MEM_verbose.out -a mem -v

    kaiju2krona -t ${nodes} -n ${names} -i ${meta.id}_kaiju_MEM_verbose.out -o ${meta.id}.krona -u

    for i in ${taxon_levels} ;
    do
      kaiju2table -t ${nodes} -n ${names} -r \$i -o ${meta.id}_kaiju_MEM.out.summary_\$i ${meta.id}_kaiju_MEM_verbose.out
    done

    grep -v U ${meta.id}_kaiju_MEM_verbose.out | gzip >  ${meta.id}_kaiju.out.gz

    rm ${meta.id}_kaiju_MEM_verbose.out

    echo \$(kaiju -h 2>&1) > v_kaiju.txt
    """
}

process KAIJU_TO_KRONA {

  publishDir "${params.outdir}/01_clean_qc/01_3_taxonomic_affiliation_reads", mode: 'copy', pattern: 'kaiju.krona.html'

  input:
    path krona_tab_files

  output:
    path "kaiju.krona.html"
    path "v_kronatools.txt", emit: v_kronatools


  script:
    """
    ktImportText -o kaiju.krona.html *.krona

    ktImportText &> v_kronatools.txt
    """
}

process PLOT_KAIJU_STAT {

  publishDir "${params.outdir}/01_clean_qc/01_3_taxonomic_affiliation_reads", mode: 'copy'

  input:
    path kaiju_results

  output:
    path "match_length_kaiju_distribution.html"

  script:
    """
    plot_kaiju_stat.py $kaiju_results -o 'match_length_kaiju_distribution.html' -v
    """
}

process MERGE_KAIJU {
  publishDir "${params.outdir}/01_clean_qc/01_3_taxonomic_affiliation_reads", mode: 'copy'

  input:
    path k_species
    path k_genus
    path k_family
    path k_class
    path k_order
    path k_phylum

  output:
    path "taxo_affi_reads_*.tsv", emit: tsv

  script:
    """
    echo '${k_species}' > species.txt
    echo '${k_genus}' > genus.txt
    echo '${k_family}' > family.txt
    echo '${k_class}' > class.txt
    echo '${k_order}' > order.txt
    echo '${k_phylum}' > phylum.txt
    for i in ${taxon_levels} ;
    do
    merge_kaiju_results.py -f \$i'.txt' -o taxo_affi_reads_\$i'.tsv'
    done
    """
}

workflow KAIJU_AND_MERGE {
  take:
    cleaned_reads
    database

  main:
    KAIJU (
      cleaned_reads,
      database
    )

    PLOT_KAIJU_STAT (
      KAIJU.out.kaiju_result.collect()
    )

    PLOT_KRONA (
      KAIJU.out.krona_tab_file.collect(), "01_clean_qc/01_3_taxonomic_affiliation_reads/", "kaiju.krona.html"
    )

    MERGE_KAIJU (
      KAIJU.out.k_species.collect(),
      KAIJU.out.k_genus.collect(),
      KAIJU.out.k_family.collect(),
      KAIJU.out.k_class.collect(),
      KAIJU.out.k_order.collect(),
      KAIJU.out.k_phylum.collect()
    )

  emit:
    report = KAIJU.out.k_all
    v_kaiju = KAIJU.out.v_kaiju
    v_kronatools = PLOT_KRONA.out.v_kronatools
}