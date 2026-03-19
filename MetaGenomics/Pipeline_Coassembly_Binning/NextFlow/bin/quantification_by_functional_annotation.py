#!/usr/bin/env python

"""----------------------------------------------------------------------------
  Script Name: quantification_by_functional_annotation.py
  Description: Aggregates per-sample read counts by functional annotation 
               categories from a gene-level abundance table. Input is a TSV 
               file with genes as rows, per-sample read count columns, and 
               functional annotation columns (GOs, KEGG_ko, KEGG_Pathway, 
               KEGG_Module, PFAMs) containing comma-separated values. For each 
               annotation type, the script: (1) explodes rows with multiple 
               annotations into separate rows (one per annotation); (2) groups 
               by annotation identifier and sums read counts across all genes 
               sharing that annotation; (3) outputs a TSV file with annotations 
               as rows, per-sample summed counts as columns, and a 'sum' column 
               with total abundance. For GO annotations only, adds a column 
               with Amigo database URLs (http://amigo.geneontology.org/amigo/term/), 
               replacing invalid entries ('-', 'nan') with '/'. Output files 
               are named <annotation>_abundance.tsv (GOs_abundance.tsv, 
               KEGG_ko_abundance.tsv, etc.) in the current directory.
  Input files: 
    - TSV file with gene-level abundances: index column = gene ID, columns = 
      per-sample read counts + functional annotation columns (GOs, KEGG_ko, 
      KEGG_Pathway, KEGG_Module, PFAMs) with comma-separated values
  Output files (in current directory):
    - GOs_abundance.tsv: GO term abundances with Amigo URL column
    - KEGG_ko_abundance.tsv: KEGG Orthology term abundances
    - KEGG_Pathway_abundance.tsv: KEGG pathway abundances
    - KEGG_Module_abundance.tsv: KEGG module abundances
    - PFAM_abundance.tsv: PFAM domain abundances
    (All outputs: annotation as index, per-sample counts, 'sum' column)
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


# Manage parameters.
parser = argparse.ArgumentParser(description = 'Aggregate per-sample read counts by functional annotation categories.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input_file', required = True, \
help = 'Name of tsv file with counts of reads and functional annotation by gene.')

parser.add_argument('-v', '--version', action = 'version', \
version = __version__)

args = parser.parse_args()

# Function of deduplication of rows when there are more than one functional annotation by column.
# From https://webdevdesigner.com/q/split-explode-pandas-dataframe-string-entry-to-separate-rows-10694/  
def splitDataFrameList(df,target_column,separator):
    ''' df = dataframe to split,
    target_column = the column containing the values to split
    separator = the symbol used to perform the split

    returns: a dataframe with each entry for the target column separated, with each element moved into a new row. 
    The values in the other columns are duplicated across the newly divided rows.
    '''
    def splitListToRows(row,row_accumulator,target_column,separator):
        split_row = str(row[target_column]).split(separator)
        for s in split_row:
            new_row = row.to_dict()
            new_row[target_column] = s
            row_accumulator.append(new_row)
    new_rows = []
    df.apply(splitListToRows,axis=1,args = (new_rows,target_column,separator))
    new_df = pd.DataFrame(new_rows)
    return new_df

# Read input file.
data = pd.read_csv(args.input_file, sep="\t", index_col=0)
print(data.head())
print(data.dtypes)
# Nb of columns.
nb_cols = data.shape[1]
data.head()
# Abundance by GOs terms.
data_GOs = splitDataFrameList(data,'GOs',',')
res_GOs = data_GOs.groupby('GOs').sum().iloc[:,0:nb_cols-23]
# Create a column to sum row abundances.
print(res_GOs)
print(res_GOs.dtypes)
res_GOs['sum'] = res_GOs.sum(axis=1)
# Create a column with amigo link.
res_GOs['GOs'] = "http://amigo.geneontology.org/amigo/term/" + res_GOs.index
res_GOs['GOs'] = res_GOs['GOs'].apply(lambda x: '/' if (x=="http://amigo.geneontology.org/amigo/term/-") else x)
res_GOs['GOs'] = res_GOs['GOs'].apply(lambda x: '/' if (x=="http://amigo.geneontology.org/amigo/term/nan") else x)
res_GOs.to_csv('GOs_abundance.tsv', sep = '\t')

# Abundance by KEGG ko.
data_KEGG_ko = splitDataFrameList(data,'KEGG_ko',',')
res_KEGG_ko = data_KEGG_ko.groupby('KEGG_ko').sum().iloc[:,0:nb_cols-23]
# Create a column to sum row abundances.
res_KEGG_ko['sum'] = res_KEGG_ko.sum(axis=1)
res_KEGG_ko.to_csv('KEGG_ko_abundance.tsv', sep = '\t')

# Abundance by KEGG_Pathway.
data_KEGG_Pathway = splitDataFrameList(data,'KEGG_Pathway',',')
res_KEGG_Pathway = data_KEGG_Pathway.groupby('KEGG_Pathway').sum().iloc[:,0:nb_cols-23]
# Create a column to sum row abundances.
res_KEGG_Pathway['sum'] = res_KEGG_Pathway.sum(axis=1)
res_KEGG_Pathway.to_csv('KEGG_Pathway_abundance.tsv', sep = '\t')

# Abundance by KEGG_Module.
data_KEGG_Module = splitDataFrameList(data,'KEGG_Module',',')
res_KEGG_Module = data_KEGG_Module.groupby('KEGG_Module').sum().iloc[:,0:nb_cols-23]
# Create a column to sum row abundances.
res_KEGG_Module['sum'] = res_KEGG_Module.sum(axis=1)
res_KEGG_Module.to_csv('KEGG_Module_abundance.tsv', sep = '\t')

# Abundance by PFAM.
data_PFAM = splitDataFrameList(data,'PFAMs',',')
res_PFAM = data_PFAM.groupby('PFAMs').sum().iloc[:,0:nb_cols-23]
# Create a column to sum row abundances.
res_PFAM['sum'] = res_PFAM.sum(axis=1)
res_PFAM.to_csv('PFAM_abundance.tsv', sep = '\t')