#!/usr/bin/env python

"""----------------------------------------------------------------------------
  Script Name: quantification_by_contig_lineage.py
  Description: Aggregates per-sample contig quantification and taxonomic 
               lineage assignments into a unified matrix where each row 
               represents a unique taxonomic lineage. For each input file 
               (typically merged idxstats + .percontig.tsv output from 
               merge_contig_quantif_perlineage.py), extracts: contig names 
               (semicolon-concatenated), contig count, read count, and mean 
               depth per lineage. Performs outer joins across samples on the 
               composite key [tax_id_by_level, lineage_by_level] to retain 
               all lineages observed in any sample. Column names are suffixed 
               with sample identifiers (derived from input filename without 
               extension) to distinguish per-sample metrics. Missing values 
               (lineages absent in a sample) are filled with 0. Output is a 
               wide-format TSV table suitable for downstream comparative 
               analysis of taxonomic abundance across samples.
  Input files: 
    - Whitespace-separated list file containing paths to per-sample TSV files 
      (format: tax_id_by_level, lineage_by_level, name_contigs, nb_contigs, 
      nb_reads, depth columns; typically <sample>_quantif_percontig.tsv)
  Output files:
    - Unified TSV table: rows = unique [tax_id_by_level, lineage_by_level] 
      pairs, columns = lineage identifiers + per-sample 
      name_contigs_<sample>/nb_contigs_<sample>/nb_reads_<sample>/depth_<sample>
      (missing values replaced with 0)
  Created By:  Mathias Vandenbogaert
  Date:        2026-05-06
-------------------------------------------------------------------------------
"""

# Metadata.
__author__ = 'Mathias Vandenbogaert'
__version__ = '0.1'
__email__ = 'mathias.vandenbogaert@lizard.bio'
__status__ = 'dev'

# Status: dev.

# Modules importation.

import argparse
import re
import sys
import pandas as pd
import os
from datetime import datetime



# Manage parameters.
parser = argparse.ArgumentParser(description = 'Aggregate per-sample lineage quantifications into a unified matrix.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--list_of_input_files', required = True, \
help = 'List of input files (one for each sample).')

parser.add_argument('-o', '--output_file', required = True, \
help = 'Name of output file containing counts of contigs and reads \
in each sample for each lineage.')

parser.add_argument('-v', '--version', action = 'version', \
version = __version__)

args = parser.parse_args()

# Recovery of the list of input files.
with open(args.list_of_input_files) as finput_list:
    sample_files = finput_list.read().split()

# Merge results for all samples by lineage.
for (sample_idx,sample_path) in enumerate(sorted(sample_files)):
    print(sample_idx)
    if(sample_idx==0):
        merge  = pd.read_csv(sample_path, delimiter='\t', dtype=str)
        sample_name = os.path.splitext(sample_path)[0]
    else:
        sample_results = pd.read_csv(sample_path, delimiter='\t', dtype=str)
        merge = pd.merge(merge,sample_results,left_on=["tax_id_by_level","lineage_by_level"],right_on=["tax_id_by_level","lineage_by_level"], how='outer', suffixes=('_' + sample_name,''))
        print (merge.head())
        sample_name = os.path.splitext(sample_path)[0]
    if('consensus_tax_id' in merge.columns): merge.drop('consensus_tax_id', inplace=True, axis=1)

# Rename columns corresponding to the last sample file.
sample_name = os.path.splitext(sample_path)[0]

merge.rename(columns = {'name_contigs': 'name_contigs_' + sample_name, \
'nb_contigs': 'nb_contigs_' + sample_name,\
'nb_reads': 'nb_reads_' + sample_name,\
'depth': 'depth_' + sample_name},inplace=True)

# Fill NaN values with 0.
merge.fillna(0, inplace=True)
print("Write " + args.output_file)
# Write merge data frame in output file.
merge.to_csv(args.output_file, sep="\t", index=False)