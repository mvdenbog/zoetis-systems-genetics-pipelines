#!/usr/bin/env python

"""----------------------------------------------------------------------------
  Script Name: merge_abundance_and_functional_annotations.py
  Description: Merges gene abundance/quantification tables with functional 
               annotations from eggNOG-mapper and Diamond best-hit results. 
               For each global gene identifier (seed_cluster), combines: 
               (1) per-sample read counts from a quantification matrix; 
               (2) functional annotations (e.g., COG, KEGG, GO, description) 
               from eggNOG-mapper output files (skipping 4 header rows); 
               (3) Diamond alignment metadata (subject accession, description) 
               from best-bitscore m8 files, with NCBI protein URL prefix added 
               to accessions and multiple hits per gene concatenated (accessions 
               with ',', descriptions with ';'). Performs left joins on gene ID 
               to preserve all quantified genes, filling missing annotation fields 
               with '-'. Outputs a unified TSV table suitable for downstream 
               functional enrichment or expression analysis.
  Input files: 
    - Counts table: TSV with header, 'seed_cluster' column + per-sample count columns
    - Annotation file list: whitespace-separated list of eggNOG-mapper TSV outputs 
      (format: standard eggNOG-mapper with 4 header rows to skip, '#query' column)
    - Diamond hits file list: whitespace-separated list of Diamond m8 best-hit files 
      (15 columns: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, 
      sstart, send, evalue, bitscore, qlen, slen, stitle; no header)
  Output files:
    - Unified TSV table: seed_cluster, per-sample counts, eggNOG annotations, 
      Diamond sseqid (NCBI URLs) and stitle (concatenated descriptions); 
      missing values replaced with '-'
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
import pandas as pd


# Manage parameters.
parser = argparse.ArgumentParser(description='Merge gene abundance tables with functional annotations and Diamond best hits.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-t', '--table_of_abundances', required=True,
                    help='Table containing counts \
for each global gene id in each sample.')

parser.add_argument('-f', '--list_of_file_annotations', required=True,
                    help='List of files storing functional annotation for each gene per sample.')

parser.add_argument('-d', '--list_of_file_diamond', required=True,
                    help='List of files storing diamond results with best bitscore \
for each gene per sample.')

parser.add_argument('-o', '--output_file', required=True,
                    help='Name of output file containing counts \
for each global gene id and its functional annotation.')

parser.add_argument('-v', '--version', action='version',
                    version=__version__)

args = parser.parse_args()

counts_file = pd.read_csv(args.table_of_abundances, header=0, sep='\t')

# Recovery of the list of annotations files.
with open(args.list_of_file_annotations) as fannot_list:
    files_of_annotations = fannot_list.read().split()

# Recovery of the list of diamond best hits files.
with open(args.list_of_file_diamond) as fdiamond_list:
    diamond_files = fdiamond_list.read().split()

# Creates a new empty dataframe for annotations.
concat_eggnog_mapper_files = pd.DataFrame()

# Concatenate annotation files.
for (annotations_idx, annotations_path) in enumerate(sorted(files_of_annotations)):
    eggnog_mapper_file = pd.read_csv(
        annotations_path, delimiter='\t', decimal='.', skiprows=4)
    concat_eggnog_mapper_files = pd.concat(
        [concat_eggnog_mapper_files, eggnog_mapper_file])

# Creates a new empty dataframe for diamond results.
concat_diamond_files = pd.DataFrame()

# Concatenate diamond files.
for (diamond_idx, diamond_path) in enumerate(sorted(diamond_files)):
    diamond_columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                       "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "stitle"]
    diamond_file = pd.read_csv(
        diamond_path, delimiter='\t', decimal='.', header=None, names=diamond_columns)
    diamond_file.loc[:, "sseqid"] = 'https://www.ncbi.nlm.nih.gov/protein/  ' + \
        diamond_file.loc[:, "sseqid"]
    group_diamond_file = diamond_file.groupby("qseqid")\
        .agg({"stitle": ';'.join, "sseqid": ','.join})\
        .reset_index()\
        .reindex(columns=diamond_file.columns)
    res_diamond_file = group_diamond_file.loc[:, [
        "qseqid", "sseqid", "stitle"]]
    concat_diamond_files = pd.concat([concat_diamond_files, res_diamond_file])

# Merge counts, annotation and diamond results.
merge_annot = pd.merge(counts_file, concat_eggnog_mapper_files,
                       left_on="seed_cluster", right_on='#query', how='left')
merge = pd.merge(merge_annot, concat_diamond_files,
                 left_on="seed_cluster", right_on="qseqid", how='left')
merge.drop('#query', inplace=True, axis=1)
merge.drop("qseqid", inplace=True, axis=1)
merge_no_nan = merge.fillna("-")
# Write merge data frame in output file.
merge_no_nan.to_csv(args.output_file, sep="\t", index=False)