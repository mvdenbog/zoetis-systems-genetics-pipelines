#!/usr/bin/env python

"""----------------------------------------------------------------------------
  Script Name: merge_kaiju_results.py
  Description: Merges Kaiju taxonomic classification summary files from multiple 
               samples into a unified abundance matrix. For a specified taxonomic 
               rank (e.g., species, genus, family), combines per-sample read 
               counts and percentages for each taxon. Input files are Kaiju 
               verbose output summaries (typically *_kaiju_MEM_verbose.out) 
               containing columns: file, taxon_id, taxon_name, reads, percent. 
               Sample names are automatically extracted by removing the 
               '_kaiju_MEM_verbose.out' suffix from the first column value. 
               Performs outer joins on taxon_name to retain all taxa observed 
               in any sample, aligning taxon_id values across files. Output is 
               a wide-format TSV table with: (1) taxon_name and taxon_id columns; 
               (2) per-sample 'reads_<sample>' and 'percent_<sample>' columns. 
               Missing values (taxa absent in a sample) are filled with 0. 
               Suitable for downstream comparative analysis or visualization 
               of taxonomic composition across samples.
  Input files: 
    - Whitespace-separated list file containing paths to Kaiju summary TSVs 
      (one per sample, same taxonomic rank; format: file, taxon_id, taxon_name, 
      reads, percent columns)
  Output files:
    - Merged TSV table: rows = taxa, columns = taxon_id + per-sample reads/percent
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

# Modules importation.

import argparse
import re
import sys
import pandas as pd
from datetime import datetime

# Manage parameters.
def parse_args():
    parser = argparse.ArgumentParser(description = 'Merge Kaiju taxonomic summary files across samples into a unified abundance matrix.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', '--list_of_kaiju_files', required = True, \
    help = 'List of kaiju summary files by the level \
    of taxonomy of interest .')

    parser.add_argument('-o', '--output_file', required = True, \
    help = 'Name of output file containing counts \
    and percentage of reads in each sample for each element \
    of the level of taxonomy of interest.')

    parser.add_argument('-v', '--version', action = 'version', \
    version = __version__)

    args = parser.parse_args()

    return args


def retrieve_annotations_files(list_kaiju_files):
    '''
    Recovery of the list of annotations files.
    '''
    with open(list_kaiju_files) as fkaiju_list:
        kaiju_files = fkaiju_list.read().split()
    return kaiju_files


def merge_kaiju_files(kaiju_files):
    '''
    Merge kaiju results for all samples.
    '''
    for (kaiju_idx,kaiju_path) in enumerate(sorted(kaiju_files)):
        if(kaiju_idx==0):
            merge  = pd.read_csv(kaiju_path, delimiter='\t', dtype=str)
        else:
            if(kaiju_idx==1):
                sample_name = merge.iloc[0,0].split('_kaiju_MEM_verbose.out')
                merge.drop('file', inplace=True, axis=1)
            else:
                sample_name = kaiju_results.iloc[0,0].split('_kaiju_MEM_verbose.out')
                merge.drop('file', inplace=True, axis=1)
            kaiju_results = pd.read_csv(kaiju_path, delimiter='\t', dtype=str)
            merge = pd.merge(merge,kaiju_results,left_on="taxon_name",\
            right_on='taxon_name', how='outer', suffixes=('_'+sample_name[0],''))
            merge['taxon_id'] = merge['taxon_id'].fillna(merge['taxon_id_' + sample_name[0]])
            merge.drop(['taxon_id_' + sample_name[0]], inplace=True, axis=1)

    # Rename columns corresponding to the last kaiju file (only if number of files > 1)
    if(kaiju_idx>0):
        sample_name = kaiju_results.iloc[0,0].split('_kaiju_MEM_verbose.out')
    else:
        sample_name = merge.iloc[0,0].split('_kaiju_MEM_verbose.out')

    merge.rename(columns = {'percent': 'percent_' + sample_name[0], \
    'reads': 'reads_' + sample_name[0]},inplace=True)

    merge.drop('file', inplace=True, axis=1)
    merge = merge[['taxon_name', 'taxon_id']  + [col for col in merge if (col != 'taxon_name' and col != 'taxon_id')]]
    # Fill the NaN by 0.
    merge.fillna(0, inplace=True)

    return merge


def main():
    args = parse_args()
    kaiju_files = retrieve_annotations_files(args.list_of_kaiju_files)
    merge = merge_kaiju_files(kaiju_files)
    merge.to_csv(args.output_file, sep="\t", index=False)


if __name__ == "__main__":
    main()