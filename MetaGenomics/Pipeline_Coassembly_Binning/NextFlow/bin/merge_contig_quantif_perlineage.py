#!/usr/bin/env python

"""----------------------------------------------------------------------------
  Script Name: merge_contig_quantif_perlineage.py
  Description: Merges per-contig read quantification (from samtools coverage) 
               with taxonomic lineage assignments (from aln_to_tax_affi.py 
               .percontig.tsv output) into aggregated abundance tables per 
               taxonomic rank. For each contig, combines: (1) mapping metrics 
               (numreads, meandepth, contig name); (2) taxonomic classification 
               (consensus lineage string, consensus taxid, taxid-by-level string). 
               Aggregates metrics by full lineage and by individual taxonomic 
               ranks (kingdom to species), computing: sum of reads, mean depth, 
               count of contigs, and semicolon-concatenated contig names. 
               Handles unmapped/unaffiliated contigs by filling NaN values with 
               placeholder values ('No affiliation', taxid=1). Outputs: (1) a 
               main TSV with lineage-level aggregations; (2) per-rank TSV files 
               for kingdom/phylum/class/order/family/genus/species; (3) Krona-
               compatible TSV files for interactive hierarchical visualization 
               of read-count and depth-based abundances. Sample name is used 
               to prefix output files and Krona outputs.
  Input files: 
    - Samtools coverage TSV: columns #rname, numreads, meandepth, endpos, etc.
    - Taxonomic affiliation TSV (.percontig.tsv): columns #contig, 
      consensus_tax_id, consensus_lineage, tax_id_by_level (semicolon-delimited)
  Output files:
    - <output_name>_quantif_percontig.tsv: aggregated metrics per full lineage
    - <output_name>_quantif_percontig_by_<rank>.tsv: aggregated metrics per 
      taxonomic rank (7 files: kingdom through species)
    - krona_reads_count_abundance/<name>.krona: Krona input TSV for read-count 
      visualization (columns: nb_reads + 7 lineage columns)
    - krona_mean_depth_abundance/<name>.krona: Krona input TSV for depth-based 
      visualization (columns: depth + 7 lineage columns)
  Created By:  Mathias Vandenbogaert
  Date:        2026-05-06
-------------------------------------------------------------------------------
"""

# Metadata.
__author__ = 'Mathias Vandenbogaert'
__version__ = '0.1'
__email__ = 'mathias.vandenbogaert@lizard.bio'
__status__ = 'dev'


from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import pandas as pd
import logging
import os

def parse_arguments():
    # Manage parameters.
    parser = ArgumentParser(description = 'Merge contig quantification and taxonomic lineage into aggregated abundance matrices per sample.',
    formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('-s', '--sam_coverage', required = True, \
    help = 'depth per contigs from samtools coverage tool.')

    parser.add_argument('-c', '--contig_tax_affi', required = True, \
    help = '.percontig.tsv file.')

    parser.add_argument("-o", '--output_name', required = True, \
    help = 'Sample name to use in output files names. The output files contains counts of contigs and reads \
    for each lineage as well as krona files.')

    parser.add_argument('-v', '--version', action = 'version', \
    version = __version__)

    parser.add_argument("--verbose", help="increase output verbosity",
                        action="store_true")

    args = parser.parse_args()
    return args

def generate_krona_directories(path, paramater,name, df_reads_count, ranks_lineage):
    '''
    path: Directory output name
    paramater: paramater analyzed (either nb_reads or depth
    '''
    read_abd_dir = path
    os.makedirs(read_abd_dir, exist_ok=True)
    outfile = os.path.join(read_abd_dir, f'{name}.krona')
    df_reads_count[[paramater ] + ranks_lineage].to_csv(outfile, index=False, header=False, sep="\t") 

def main():

    args = parse_arguments()

    if args.verbose:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG)
        logging.info('Mode verbose ON')

    else:
        logging.basicConfig(format="%(levelname)s: %(message)s")

    sam_coverage_file = args.sam_coverage
    contig_taxaffi_file = args.contig_tax_affi
    name =  args.output_name
    output_name = f'{name}_quantif_percontig' 

    ranks = ["kingdom", "phylum", "class", "order",
            "family", "genus", "species"]

    logging.info("Read and merge tables")
    cov_df = pd.read_csv(sam_coverage_file, delimiter='\t')

    contig_taxaffi_df = pd.read_csv(contig_taxaffi_file, delimiter='\t', dtype=str)

    depth_tax_contig_df = pd.merge(cov_df,contig_taxaffi_df,left_on='#rname',right_on='#contig', how='outer')

    # Fill NaN values to keep unmapped contigs.
    depth_tax_contig_df['consensus_lineage'] = depth_tax_contig_df['consensus_lineage'].fillna('No affiliation')
    depth_tax_contig_df['tax_id_by_level'] = depth_tax_contig_df['tax_id_by_level'].fillna(1)
    depth_tax_contig_df['consensus_tax_id'] = depth_tax_contig_df['consensus_tax_id'].fillna(1)


    logging.info("group by lineage")
    groupby_cols = ['consensus_lineage','consensus_tax_id', 'tax_id_by_level']
    depth_lineage_df = depth_tax_contig_df.groupby(groupby_cols).agg({
                                                '#rname' : [';'.join, 'count'],
                                                'numreads': 'sum',
                                                'meandepth': 'mean'}).reset_index()

    depth_lineage_df.columns=['lineage_by_level', 'consensus_tax_id', 'tax_id_by_level',
                              'name_contigs', 'nb_contigs', 'nb_reads', 'depth']

    logging.info(f"Write out {output_name}.tsv")
    depth_lineage_df.to_csv(f"{output_name}.tsv", sep="\t", index=False)


    # split lineage
    ranks_taxid = [f"{r}_taxid" for r in ranks]
    ranks_lineage = [f"{r}_lineage" for r in ranks]

    try:
        depth_lineage_df[ranks_taxid] = depth_lineage_df['tax_id_by_level'].str.split(pat=";",expand=True)
        depth_lineage_df[ranks_lineage] = depth_lineage_df["lineage_by_level"].str.split(pat=";",expand=True)

    except (ValueError,AttributeError):
        # Manage case when lineage_by_level is only equal to " Unable to find taxonomy consensus" or "No affiliation"
        df_noaffi = pd.DataFrame("no_affi", index=range(len(depth_lineage_df)), columns=ranks_taxid+ranks_lineage)
        depth_lineage_df  = pd.concat([depth_lineage_df, df_noaffi], axis=1)
    depth_lineage_df = depth_lineage_df.fillna(value='no_affi')

    # groupby each rank and write the resulting table
    levels_columns=['tax_id_by_level','lineage_by_level','name_contigs','nb_contigs', 'nb_reads', 'depth']

    logging.info("group by rank")
    for rank in ranks:
        depth_rank_lineage_df = depth_lineage_df.groupby([f'{rank}_taxid',f'{rank}_lineage']).agg({
                                                      'name_contigs' : [';'.join],
                                                      'nb_contigs' : 'sum',
                                                      'nb_reads' : 'sum',
                                                      'depth': 'mean'}).reset_index()

        depth_rank_lineage_df.columns=levels_columns
        depth_rank_lineage_df['rank'] = rank
        logging.info(f"Write out {output_name}_by_{rank}.tsv")
        depth_rank_lineage_df.to_csv(f"{output_name}_by_{rank}.tsv", sep="\t", index=False)

    logging.info('Writting krona tsv files')
    # create krona tsv file used to generate krona plot
    df_reads_count = depth_lineage_df[["nb_reads", "depth"] + ranks_lineage]
    df_reads_count = df_reads_count.replace(' None', "") 

    # number of reads abundance
    generate_krona_directories("krona_reads_count_abundance", "nb_reads", name, df_reads_count, ranks_lineage)
    # mean depth abundance
    generate_krona_directories("krona_mean_depth_abundance", "depth", name, df_reads_count, ranks_lineage)

if __name__ == '__main__':
    main()