#!/usr/bin/env python3

"""----------------------------------------------------------------------------
  Script Name: rename_contigs.py
  Description: Renames contig identifiers in a FASTA assembly file to a 
               standardized format: <sample_name>_c<contig_number>, where 
               contig_number is a sequential integer (1-based) assigned in 
               the order contigs appear in the input file. The original 
               contig name is preserved in the FASTA header comment field 
               (after the new ID). Additionally generates a two-column TSV 
               mapping file recording the original and new contig names for 
               traceability. Output FASTA uses standard format with sequence 
               lines unchanged. Supports any FASTA input readable by pyfastx. 
               Useful for standardizing contig naming across multiple samples 
               prior to binning, annotation, or comparative analysis.
  Input files: 
    - Input FASTA file: assembly contigs with arbitrary identifiers
  Output files:
    - Renamed FASTA file: contigs with IDs formatted as <sample>_c<N>, 
      original ID preserved in header comment
    - Mapping TSV file (default: original_to_new_contig_name.tsv): two 
      columns (original_name, new_name) for reverse lookup
  Created By:  Mathias Vandenbogaert
  Date:        2026-05-06
-------------------------------------------------------------------------------
"""

# Metadata
__author__ = 'Mathias Vandenbogaert'
__version__ = '0.1'
__email__ = 'mathias.vandenbogaert@lizard.bio'
__status__ = 'dev'


from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
import logging
import gzip
import csv
from collections import defaultdict
import pyfastx
from Bio.Seq import Seq


def parse_arguments():
    """Parse script arguments."""
    parser = ArgumentParser(description="Rename contig identifiers to standardized sample_cN format.",
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('-s', '--sample', help='Sample name used to rename contigs.', required=True)
    
    parser.add_argument('-i', '--fna_file', help='Original fasta file of contigs.', required=True)

    parser.add_argument('-o', '--out_fna', help='Output fasta file with renamed contigs.', required=True)
    
    parser.add_argument('-t', '--contig_names_table', help='Tabular table with 2 fields : orignal and new name.', default="original_to_new_contig_name.tsv")

    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")

    args = parser.parse_args()
    return args



def main():

    args = parse_arguments()

    if args.verbose:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG)
        logging.info('Mode verbose ON')

    else:
        logging.basicConfig(format="%(levelname)s: %(message)s")

    sample = args.sample    
    fna_file = args.fna_file
    out_fna = args.out_fna


    logging.info(f'Writting renamed fasta file in {out_fna}')
    with open(out_fna, "w") as fh_fna:
        for i, (name, seq) in enumerate(pyfastx.Fasta(fna_file, build_index=False)):
            fh_fna.write(f'>{sample}_c{i+1} {name}\n{seq}\n')

    sample = args.sample    
    fna_file = args.fna_file
    out_fna = args.out_fna
    old2new_contig_name = args.contig_names_table


    logging.info(f'Writting renamed fasta file in {out_fna}')
    with open(out_fna, "w") as fh_fna, open(old2new_contig_name, 'w') as mout:
        for i, (name, seq) in enumerate(pyfastx.Fasta(fna_file, build_index=False)):
            fh_fna.write(f'>{sample}_c{i+1} {name}\n{seq}\n')
            mout.write(f'{name}\t{sample}_c{i+1}\n')

if __name__ == '__main__':
    main()