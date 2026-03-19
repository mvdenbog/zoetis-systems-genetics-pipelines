#!/usr/bin/env python

"""--------------------------------------------------------------------
  Script Name: merge_kaiju_results.py
  Description: join kaiju results by taxon level for all samples.
  Input files: All kaiju summary files by level of taxonomy of interest.
-----------------------------------------------------------------------
"""

# Metadata.
__author__ = 'Mathias Vandenbogaert'
__version__ = '0.1'
__status__ = 'dev'

# Modules importation.

import argparse
import re
import sys
import pandas as pd
from datetime import datetime

# Manage parameters.
def parse_args():
    parser = argparse.ArgumentParser(description = 'Script which join \
    kaiju results by level of taxonomy of interest for all samples.')

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
