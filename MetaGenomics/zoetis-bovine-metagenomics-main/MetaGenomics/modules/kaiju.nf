include {PLOT_KRONA} from  './krona'

taxon_levels = "phylum class order family genus species"

process KAIJU {
disk 150.GB
  tag "${meta.id}"
  label "SPADES"

  publishDir "${params.outdir}/${meta.id}/Kaiju", mode: 'copy', pattern: '*_kaiju.out.gz'
  publishDir "${params.outdir}/${meta.id}/Kaiju", mode: 'copy', pattern: '*.fq'
  publishDir "${params.outdir}/${meta.id}/Kaiju", mode: 'copy', pattern: '*.krona'
  publishDir "${params.outdir}/${meta.id}/Kaiju", mode: 'copy', pattern: '*_kaiju_krona.html'

  input:
    tuple val(meta), path(reads)
    tuple path(nodes), path(fmi), path(names)

  output:
    path  "${meta.id}_kaiju.out.gz",  emit: kaiju_result
    path "${meta.id}.krona", emit: krona_tab_file
    path "${meta.id}_concat.fastq"
    path "*.list"
    path "*.summary_*", emit: k_all
    path "*.summary_species", emit: k_species
    path "*.summary_genus", emit: k_genus
    path "*.summary_family", emit: k_family
    path "*.summary_class", emit: k_class
    path "*.summary_order", emit: k_order
    path "*.summary_phylum", emit: k_phylum
    path "v_kaiju.txt", emit: v_kaiju
    path "*_kaiju_labels.out", emit: kaiju_labels
    tuple val("${meta.id}"), path("*.fq"), emit: kaiju_fq_species
    path "*_kaiju_krona.html"
  script:
    """

zcat ${reads[0]} ${reads[1]} > ${meta.id}_concat.fastq


kaiju -z ${task.cpus} -t ${nodes} -f ${fmi} -i ${meta.id}_concat.fastq  -o ${meta.id}_kaiju.out -v -a greedy -m 1 -e 5
kaiju-addTaxonNames -i ${meta.id}_kaiju.out -o ${meta.id}_kaiju_labels.out -t ${nodes} -n ${names} -u -p
kaiju2krona -i ${meta.id}_kaiju.out -t ${nodes} -n ${names} -o ${meta.id}.krona
ktImportText -o  ${meta.id}_kaiju_krona.html ${meta.id}.krona

    for i in ${taxon_levels} ;
    do
      kaiju2table -t ${nodes} -n ${names} -r \$i -o ${meta.id}_kaiju.out.summary_\$i ${meta.id}_kaiju.out
    done

    grep -v U ${meta.id}_kaiju.out | gzip >  ${meta.id}_kaiju.out.gz


# seqkit  fq2fa ${meta.id}_concat.fastq -o out.fasta

for species in Virus Bacteria Fungi; do
grep -i \${species} ${meta.id}_kaiju_labels.out | cut -f 2 | sort | uniq > \${species}.list
get_fastq_from_prefix_seq.pl \${species}.list ${meta.id}_concat.fastq ${meta.id}.\${species}.fq
done

grep '^C' ${meta.id}_kaiju_labels.out | grep -v -i virus | grep -v -i bacteria | grep -v -i fungi | grep -v -i eukaryota | cut -f 2 | sort | uniq > ${meta.id}.Other.list
get_fastq_from_prefix_seq.pl ${meta.id}.Other.list ${meta.id}_concat.fastq ${meta.id}.Other.fq



grep '^U' ${meta.id}_kaiju.out | cut -f 2 | sort | uniq > ${meta.id}.unclassifieds.list
get_fastq_from_prefix_seq.pl ${meta.id}.unclassifieds.list  ${meta.id}_concat.fastq ${meta.id}.unclassifieds.fq

    echo \$(kaiju -h 2>&1) > v_kaiju.txt
    """
}

process KAIJU_TO_KRONA {
label "GENERAL"

  publishDir "${params.outdir}/Kaiju/", mode: 'copy', pattern: 'kaiju.krona.html'

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
label "GENERAL"

  publishDir "${params.outdir}/Kaiju", mode: 'copy'

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
label "GENERAL"
  publishDir "${params.outdir}/Kaiju/", mode: 'copy'

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

//    PLOT_KAIJU_STAT (
//      KAIJU.out.kaiju_result.collect()
//    )

    PLOT_KRONA (
      KAIJU.out.krona_tab_file.collect(), "Kaiju/", "kaiju.krona.html"
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
    kaiju_labels = KAIJU.out.kaiju_labels
    kaiju_fq_species = KAIJU.out.kaiju_fq_species
}
