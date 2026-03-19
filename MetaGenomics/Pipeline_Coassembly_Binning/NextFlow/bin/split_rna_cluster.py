#!/usr/bin/env python

"""----------------------------------------------------------------------------
  Script Name: split_rna_clusters.py
  Description: Separates a featureCounts gene-level count table into three 
               distinct files based on RNA feature type: (1) ribosomal RNA 
               (rRNA) genes; (2) transfer RNA (tRNA) genes; (3) protein-coding 
               sequences (CDS). Input is a featureCounts output TSV where the 
               first column contains gene identifiers. Genes with 'rRNA' or 
               'tRNA' in their identifier are classified as RNA; all others 
               are treated as CDS. The rRNA/tRNA subset is further split by 
               detecting 'tRNA' to separate tRNA from rRNA. All three output 
               files preserve the original two-line header from featureCounts. 
               A post-processing step removes a duplicated first data row that 
               can occur in the CDS output due to pandas CSV writing behavior. 
               Output filenames are prefixed: cds.<input>, rrna.<input>, 
               trna.<input>. Useful for separating RNA and coding gene counts 
               prior to differential expression or functional analysis.
  Input files: 
    - featureCounts TSV: gene count table with gene identifiers in first column
  Output files (in current directory):
    - cds.<input_filename>: protein-coding gene counts
    - rrna.<input_filename>: ribosomal RNA gene counts
    - trna.<input_filename>: transfer RNA gene counts
    (All outputs preserve original 2-line header, tab-separated format)
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
from argparse import ArgumentParser
import pandas as pd

# Manage parameters.
def parse_arguments():
    """Parse script arguments."""
    parser = ArgumentParser(description='Separate the rrna, trna and cds lines into 3 different files')

    parser.add_argument('-f', '--file_of_counts', required = True,
                        help = 'Table of counts for each gene from featureCounts.')
    
    args = parser.parse_args()
    return args

### Functions
def split_rna_cluster(file_of_counts):
    # Read count table
    df = pd.read_csv(file_of_counts)
    # Find rna lines
    condition = df.iloc[:, 0].astype(str).str.contains('tRNA|rRNA')
    cds = df[-condition]
    rna = df[condition]
    
    conditionT = rna.iloc[:, 0].astype(str).str.contains('tRNA')
    rrna = rna[-conditionT]
    trna = rna[conditionT]

    with open(file_of_counts, 'r') as f:
        headers = [next(f) for _ in range(2)]

    # Output path + file name
    cds_name = "cds." + file_of_counts
    rrna_name = "rrna." + file_of_counts
    trna_name = "trna." + file_of_counts

    with open(cds_name, 'w') as f_cds, open(rrna_name, 'w') as f_rrna, open(trna_name, 'w') as f_trna:
        f_cds.writelines(headers)
        f_rrna.writelines(headers)
        f_trna.writelines(headers)

    # Write splitted tables
    cds.to_csv(cds_name, mode='a', index=False, header=False)
    rrna.to_csv(rrna_name, mode='a', index=False, header=False)
    trna.to_csv(trna_name, mode='a', index=False, header=False)

    # Remove duplicated first row in cds table
    with open(cds_name, 'r') as f:
        lines = f.readlines()
    
    if len(lines) > 1:
        del lines[1]

    with open(cds_name, 'w') as f:
        f.writelines(lines)

def main():
    args = parse_arguments()
    split_rna_cluster(args.file_of_counts)

if __name__ == "__main__":
    main()