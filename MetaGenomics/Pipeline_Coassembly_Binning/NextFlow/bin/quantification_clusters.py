#!/usr/bin/env python

"""----------------------------------------------------------------------------
  Script Name: quantification_clusters.py
  Description: Aggregates gene-level read counts to global cluster-level 
               abundances by resolving a two-tier clustering hierarchy. 
               Input consists of: (1) a table mapping global cluster IDs 
               (seed clusters from final CD-HIT) to intermediate cluster IDs 
               (from per-sample CD-HIT runs); (2) a list of files containing 
               intermediate cluster ID to gene ID mappings per sample; 
               (3) a list of featureCounts output files with per-gene read 
               counts per sample. The script builds a global cluster -> gene 
               correspondence table by chaining the two mapping levels, then 
               sums read counts from featureCounts files into the appropriate 
               global cluster for each sample. Outputs: (1) a TSV file mapping 
               each gene to its global cluster ID (seed_cluster); (2) a TSV 
               matrix with global clusters as rows, samples as columns, and 
               aggregated read counts as values. Gene IDs are extracted from 
               featureCounts output by removing the '_gene' suffix.
  Input files: 
    - Correspondence table (TSV): global_cluster_id <tab> intermediate_cluster_id(s)
    - Cluster-gene list file: whitespace-separated paths to files containing 
      intermediate_cluster_id <tab> gene_id mappings (one file per sample)
    - Counts list file: whitespace-separated paths to featureCounts output 
      files (.featureCounts.count format with Geneid column and count in 7th field)
  Output files:
    - Global cluster <-> gene mapping TSV: columns seed_cluster, id_gene
    - Global cluster abundance matrix TSV: columns seed_cluster + per-sample 
      count columns (sample names derived from input file paths)
  Created By:  Mathias Vandenbogaert
  Date:        2026-05-06
-------------------------------------------------------------------------------
"""

# Metadata.
__author__ = 'Mathias Vandenbogaert'
__version__ = '0.1'
__email__ = 'mathias.vandenbogaert@lizard.bio'
__status__ = 'dev'

# Modules importation.
import argparse
import logging
from datetime import datetime

### Functions
def parse_arguments():
    '''
    Parse parameters.
    '''
    parser = argparse.ArgumentParser(description = 'Aggregate gene-level read counts to global cluster-level abundances via two-tier clustering hierarchy.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-t', '--table_of_corespondences', required = True, 
    help = 'Correspondence table between global cluster \
    id and intermediate cluster id.')

    parser.add_argument('-l', '--list_of_file_clusters', required = True, 
    help = 'List of files containing correspondence tables between \
    cluster intermediate cluster id and gene id per sample.')

    parser.add_argument('-c', '--list_of_file_counts', required = True, 
    help = 'List of files storing read counts for each gene per sample.')

    parser.add_argument('-oc', '--output_counts', required = True, 
    help = 'Name of output file containing counts \
    for each global cluster id and each sample.')

    parser.add_argument('-oid', '--output_id', required = True, 
    help = 'Name of output file containing correspondence table \
    between global cluster id and gene id.')

    parser.add_argument("--verbose", help="increase output verbosity",
    action="store_true")

    parser.add_argument('-v', '--version', action = 'version', \
    version = __version__)

    args = parser.parse_args()
    return args

# Recovery of the list of file names.
def processing_input_files(file_counts, correspondance_table, file_clusters):
    '''
    Recovering the list of file names.
    For all variable names:
    g_clstr: global cluster,
    int_clstr: intermediate cluster,
    gene: gene.
    
    '''
    with open(file_counts) as fcounts_list:
        files_of_counts = fcounts_list.read().split()


    d_g_clstr_id_by_int_clstr_id = {}
    d_count_by_g_clstr = {}

    with open(correspondance_table) as fp:
        for g_clstr_int_clstr_line in fp:
            g_clstr, *int_clstr = g_clstr_int_clstr_line.split()
            for clstr in int_clstr :
                d_g_clstr_id_by_int_clstr_id[clstr] = g_clstr
            d_count_by_g_clstr[g_clstr] = [0]*len(files_of_counts)
    
    d_g_clstr_id_by_gene_id = {}

    # Store into files_of_int_clstr_id_gene_id the list of sample files names
    # which contains correspondence between intermediate cluster id and gene id.
    #
    # For each line of each sample file into files_of_int_clstr_id_gene_id,
    # store the gene id (key) in the dictionnary
    # d_g_clstr_id_by_gene_id.
    # The value of d_g_clstr_id_by_gene_id is the value of
    # d_g_clstr_id_by_int_clstr_id (global cluster id).

    with open(file_clusters) as fcluster_list:
        files_of_int_clstr_id_gene_id = fcluster_list.read().split()

    for int_clstr_gene_path in files_of_int_clstr_id_gene_id:
        with open(int_clstr_gene_path) as fh:
            for file_int_clstr_gene in fh:
                line_int_clstr_gene = file_int_clstr_gene.split()
                int_clstr_id = line_int_clstr_gene[0]
                gene_id_from_clstr_gene_path = line_int_clstr_gene[1]
                if 'd_g_clstr_id_by_gene_id[gene_id_from_clstr_gene_path]' not in d_g_clstr_id_by_gene_id:
                    d_g_clstr_id_by_gene_id[gene_id_from_clstr_gene_path] \
                    = d_g_clstr_id_by_int_clstr_id[int_clstr_id]
                else:
                    d_g_clstr_id_by_gene_id[gene_id_from_clstr_gene_path]\
                    .append(d_g_clstr_id_by_int_clstr_id[int_clstr_id])

    return files_of_counts, d_count_by_g_clstr, d_g_clstr_id_by_gene_id


def linking_counts_and_clusters(files_of_counts, d_count_by_g_clstr, d_g_clstr_id_by_gene_id):
    '''
    For each count file (output of featureCounts), reading of lines one by one,
    recovery of name of gene and count number and incrementing of corresponding
    value in d_count_by_g_clstr.
    '''
    for (count_idx,counts_path) in enumerate(files_of_counts):
        with open(counts_path) as fh:
            for f_gene_counts in fh:
                if f_gene_counts.startswith('#') \
                or f_gene_counts.startswith('Geneid'):
                    continue
                line_gene_counts_split = f_gene_counts.split()
                gene_id = line_gene_counts_split[0].split("_gene")[0]
                gene_count = int(line_gene_counts_split[6])
                d_count_by_g_clstr[d_g_clstr_id_by_gene_id[gene_id]]\
                [count_idx] += gene_count
    
    return d_count_by_g_clstr

def writing_outputs( files_of_counts, d_count_by_g_clstr, d_g_clstr_id_by_gene_id, output_id, output_counts):
    '''
    # Write output file containing correspondence table
    # between global cluster id and gene id.
    '''
    with open(output_id,"w") as foutput_res_table:
        # Heading of output file: name of columns.
        foutput_res_table.write("seed_cluster" + "\t" + "id_gene" + "\n")
        # Writing seed cluster ids and genes ids for each sample contained in
        # d_g_clstr_id_by_gene_id in the output file line by line.
        for gene_id, g_clstr_id \
        in d_g_clstr_id_by_gene_id.items():
            foutput_res_table.write(g_clstr_id \
            + "\t" \
            + gene_id \
            + "\n")

    # Write output file containing global cluster id and read count for each sample.
    with open(output_counts,"w") as foutput_res_counts:
        # Heading of output file: name of columns.
        foutput_res_counts.write("seed_cluster\t" + "\t".join(files_of_counts) + "\n")
        # Writing global cluster ids and counts for each sample contained in
        # d_count_by_g_clstr in the output file line by line.
        for g_clstr, count in d_count_by_g_clstr.items():
            foutput_res_counts.write(g_clstr + "\t" \
            + "\t".join([str(i) for i in count])\
            + "\n")

def main():
    args = parse_arguments()

    if args.verbose:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG)
        logging.info('Mode verbose ON')

    else:
        logging.basicConfig(format="%(levelname)s: %(message)s")

    # Print time.
    logging.info(str(datetime.now()))

    files_of_counts, d_count_by_g_clstr, d_g_clstr_id_by_gene_id = \
    processing_input_files(args.list_of_file_counts, args.table_of_corespondences, args.list_of_file_clusters)

    d_count_by_g_clstr = linking_counts_and_clusters(files_of_counts, d_count_by_g_clstr, d_g_clstr_id_by_gene_id)

    writing_outputs( files_of_counts, d_count_by_g_clstr, d_g_clstr_id_by_gene_id, args.output_id, args.output_counts)

if __name__ == '__main__':
    main()