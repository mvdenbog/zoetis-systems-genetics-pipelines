include { ASSIGN_TAXONOMY;PLOT_TAXONOMIC_AFFILIATIONS } from '../modules/assign_taxonomy'
include { QUANTIF_AND_TAXONOMIC_TABLE_CONTIGS } from '../modules/quantif_and_taxonomic_table_contigs'
include {PLOT_KRONA as KRONA_READS_COUNT; PLOT_KRONA as KRONA_DEPTH} from  '../modules/krona'

workflow STEP_07_TAXO_AFFI {
    take:
        taxonomy
        diamond_result // channel: [ val(meta), path(diamond_file) ]
        sam_coverage // channel: [ val(meta), path(samtools coverage) ]
    main:
        if (params.coassembly){
            sam_coverage.map { meta, cov ->
                             [ meta.group, meta, cov ]}
                        .set { ch_sam_cov_tmp }

            diamond_result.map { meta, m8 ->
                               [ meta.group, m8]} 
                          .combine( ch_sam_cov_tmp, by: 0)
                          .map { group, m8, meta, cov -> 
                               [ meta, m8, cov ]}
                          .set { ch_assign_taxo_input }
        } else {
            ch_assign_taxo_input = diamond_result.join(sam_coverage, remainder: true)
        }

        ASSIGN_TAXONOMY ( taxonomy, ch_assign_taxo_input )

        QUANTIF_AND_TAXONOMIC_TABLE_CONTIGS (
            ASSIGN_TAXONOMY.out.q_all.collect(),
            ASSIGN_TAXONOMY.out.q_superkingdom.collect(),
            ASSIGN_TAXONOMY.out.q_phylum.collect(),
            ASSIGN_TAXONOMY.out.q_order.collect(),
            ASSIGN_TAXONOMY.out.q_class.collect(),
            ASSIGN_TAXONOMY.out.q_family.collect(),
            ASSIGN_TAXONOMY.out.q_genus.collect(),
            ASSIGN_TAXONOMY.out.q_species.collect()
        )

        PLOT_TAXONOMIC_AFFILIATIONS(ASSIGN_TAXONOMY.out.q_all.collect())

        KRONA_READS_COUNT(ASSIGN_TAXONOMY.out.krona_reads_count.collect(), "07_taxo_affi/07_3_plot", "krona_read_count_abundance.html")

        KRONA_DEPTH(ASSIGN_TAXONOMY.out.krona_depth.collect(), "07_taxo_affi/07_3_plot", "krona_mean_depth_abundance.html")

    emit:
        quantif_by_contig_lineage = QUANTIF_AND_TAXONOMIC_TABLE_CONTIGS.out.quantif_by_contig_lineage
        software_versions = KRONA_DEPTH.out.v_kronatools
}
