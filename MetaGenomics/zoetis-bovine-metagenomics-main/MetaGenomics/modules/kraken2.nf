process KRAKEN2 {

 publishDir "${params.outdir}/${meta.id}/Kraken2/", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path(database)

    output:
    tuple val("kraken2"), val(meta), path("results.krona"), emit: results_for_krona
    tuple val(meta), path("*kraken2_report.txt")          , emit: report
    path "versions.yml"                                   , emit: versions
    path "v_kraken2.txt", emit: v_kraken2

    script:
    def input = meta.single_end ? "\"${reads}\"" :  "--paired \"${reads[0]}\" \"${reads[1]}\""
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    kraken2 \
        --report-zero-counts \
        --threads ${task.cpus} \
        --db ${database} \
        --report ${prefix}.kraken2_report.txt \
        $input \
        > kraken2.kraken
    cat kraken2.kraken | cut -f 2,3 > results.krona

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //' | sed 's/ Copyright.*//')
    END_VERSIONS
    kraken2 --version 2>&1 > v_kraken2.txt
    """
}

process KRAKEN2ToKRONA{
publishDir "${params.outdir}/${meta.id}/Kraken2/Krona/", mode: 'copy'
input:
tuple val(meta), path(kraken2_report)
path (k2_taxonomy_db)

output:
path "${meta.id}_krona.in", emit: krona_in
path "*.html"

script:
"""
cut -f 2,5 ${kraken2_report} > ${meta.id}_krona.in
ktImportTaxonomy -k -tax ${k2_taxonomy_db} -m 1 -o ${meta.id}_kraken2_krona.html ${meta.id}_krona.in

"""
}

process KRAKEN2KRONAMERGE{
publishDir "${params.outdir}/Kraken2/KronaMerged/", mode: 'copy'
// container = 'europe-west1-docker.pkg.dev/zoetis-bovine-metagenomics/containers/kaiju_seqkit_spades:1.2'
input:
path(files)
path(k2_taxonomy_db)
output:
path "kraken2_krona.html"
script:
"""
ktImportTaxonomy -k -tax ${k2_taxonomy_db} -m 1 -o kraken2_krona.html *.in
"""
}


