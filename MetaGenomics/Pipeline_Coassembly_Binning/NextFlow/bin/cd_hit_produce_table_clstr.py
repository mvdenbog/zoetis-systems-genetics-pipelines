#!/usr/bin/env python

"""----------------------------------------------------------------------------
  Script Name: cdhit_clusters_to_table.py
  Description: Parses CD-HIT clustering output (.clstr format) and converts it 
               into a two-column TSV table mapping each cluster representative 
               sequence to all its member sequences (including itself). For each 
               CD-HIT cluster, the representative sequence (marked with '*' in 
               the .clstr file) is listed in the first column, and every sequence 
               belonging to that cluster (representative + members) is listed 
               in the second column, one row per member. This format facilitates 
               downstream analysis of cluster membership, dereplication tracking, 
               or abundance aggregation by cluster.
  Input files: 
    - CD-HIT .clstr output file (standard format with '>' headers, '*' marking 
      representative sequences, and '...>' delimiters for sequence names)
  Output files:
    - TSV table with two columns: 'representative' and 'member', where each row 
      associates one cluster member with its cluster's representative sequence
  Created By:  Mathias Vandenbogaert
  Date:        2026-05-06
-------------------------------------------------------------------------------
"""

# Metadata
__author__ = 'Mathias Vandenbogaert'
__version__ = '0.1'
__email__ = 'mathias.vandenbogaert@lizard.bio'
__status__ = 'dev'

import sys
import re
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def parse_arguments():
    """Parse script arguments."""
    parser = ArgumentParser(
        description="Convert CD-HIT .clstr clustering output to representative-member TSV table.",
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-i', '--input_file',
        required=True,
        help="CD-HIT output file (.clstr format) representing clusters."
    )
    parser.add_argument(
        '-o', '--output_file',
        required=True,
        help="Output TSV table mapping representatives to cluster members."
    )
    args = parser.parse_args()
    return args


def process(input_file, output_file):
    """
    Parse CD-HIT .clstr file and write representative-member mappings to TSV.
    
    CD-HIT .clstr format:
    - Cluster headers start with '>'
    - Each sequence line contains: index>length, >seq_name... at '*' if representative
    - Representative sequences are marked with '*' at the end of the line
    
    Output: TSV with columns [representative, member], one row per cluster member.
    """
    # Initialize tracking variables
    ref = ""
    seqs = []
    
    with open(input_file, 'r') as FH_input, open(output_file, 'w') as FH_out:
        for line in FH_input:
            if line == '':
                break
            else:
                if line[0] == ">":
                    # Write previous cluster mappings before starting new cluster
                    for seq in seqs:
                        FH_out.write(ref + "\t" + seq + "\n")
                    # Reset for new cluster
                    ref = ""
                    seqs = []
                else:
                    # Parse sequence line: extract name and check if representative
                    a, b = line.split('>', 1)
                    name = b.split("...")[0]
                    rep = (line.rstrip()[-1] == '*')
                    if rep:
                        ref = name
                        seqs.append(name)
                    else:
                        seqs.append(name)
        # Write final cluster mappings after loop ends
        for seq in seqs:
            FH_out.write(ref + "\t" + seq + "\n")


def main():
    """Main execution function."""
    args = parse_arguments()
    process(args.input_file, args.output_file)


if __name__ == '__main__':
    main()