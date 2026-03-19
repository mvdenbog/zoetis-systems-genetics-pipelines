#!/usr/bin/env python

"""----------------------------------------------------------------------------
  Script Name: write_unbinned_contigs.py
  Description: Extracts contigs from a metagenomic assembly that were not 
               assigned to any final bin during the binning process. The script 
               reads one or more bin FASTA files to collect all binned contig 
               identifiers, then iterates through the original assembly FASTA 
               and writes only those contigs not found in any bin to a separate 
               output file. This enables downstream analysis of unbinned 
               sequence content, assessment of binning completeness, or 
               recovery of potentially novel genomic elements excluded from 
               bins. Uses pyfastx for efficient FASTA parsing without full 
               indexing.
  Input files: 
    - Assembly FASTA: original metagenomic assembly containing all contigs
    - Bin FASTA files (one or more): final bin sequences from binning tools 
      (e.g., MetaBAT2, MaxBin, CONCOCT, or dereplicated sets)
  Output files:
    - unbinned_contigs.fasta (default): FASTA file containing only contigs 
      not present in any input bin file
  Created By:  Mathias Vandenbogaert
  Date:        2026-05-06
-------------------------------------------------------------------------------
"""

# Metadata
__author__ = 'Mathias Vandenbogaert'
__version__ = '0.1'
__email__ = 'mathias.vandenbogaert@lizard.bio'
__status__ = 'dev'


import argparse
import pyfastx


def retrieve_contigs_in_bins(bins_files):
    binned_contigs = set()
    for bin_fasta in bins_files:
        binned_contigs |= set(pyfastx.Fasta(bin_fasta).keys())
    return binned_contigs


def write_unbinned_contigs(binned_contigs, assembly_file, unbinned_outfile):

    with open(unbinned_outfile, 'w') as out_fh:
        for contig, seq in pyfastx.Fasta(assembly_file, build_index=False):
            if contig not in binned_contigs:
                out_fh.write(f'>{contig}\n{seq}\n')


def parse_arguments():
    # Manage parameters.
    parser = argparse.ArgumentParser(description='Extract unbinned contigs from assembly based on bin membership.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-b', '--bins', nargs='+', required = True, 
                        help = 'List of path samples bins to check')

    parser.add_argument('-a', '--assembly', required = True, 
                        help = 'Path to assembly')

    parser.add_argument('-o', '--unbinned', default='unbinned_contigs.fasta',
                        help = 'Output file with unbinned contigs file.')

    args = parser.parse_args()
    return args

def main():

    args = parse_arguments()

    binned_contigs = retrieve_contigs_in_bins( args.bins)
    write_unbinned_contigs(binned_contigs, args.assembly, args.unbinned )
    

if __name__ == '__main__':
    main()