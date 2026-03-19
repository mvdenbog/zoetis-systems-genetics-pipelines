#!/usr/bin/env python

"""----------------------------------------------------------------------------
  Script Name: Filter_contig_per_cpm.py
  Description: Filters metagenomic contigs based on CPM (Counts Per Million) 
               normalization of mapped reads. For each contig, CPM is calculated 
               as: 1e6 * (mapped_reads / total_mapped_reads_across_all_contigs). 
               Contigs with CPM >= the specified cutoff are retained. Supports 
               processing multiple samtools idxstats files (e.g., from multiple 
               samples): mapped read counts and CPM values are summed across 
               files for each contig, and a contig is kept if it passes the 
               cutoff in at least one sample (union logic). Outputs two FASTA 
               files: one containing retained contigs and one containing 
               discarded contigs. Uses pyfastx for efficient FASTA parsing.
  Input files: 
    - One or more samtools idxstats TSV files (columns: ref_name, length, 
      mapped_reads, unmapped_reads; '*' comment lines ignored)
    - Reference FASTA file containing contig sequences
  Output files:
    - Selected FASTA: contigs with CPM >= cutoff (specified via --select)
    - Discarded FASTA: contigs with CPM < cutoff (specified via --discard)
  Created By:  Mathias Vandenbogaert
  Date:        2026-05-06
-------------------------------------------------------------------------------
"""

# Metadata
__author__ = 'Mathias Vandenbogaert'
__version__ = '0.1'
__email__ = 'mathias.vandenbogaert@lizard.bio'
__status__ = 'dev'

# Status: dev

# Modules importation


from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import pandas as pd
import logging
import pyfastx

################################################
# Function
################################################

def parse_arguments():
    """Parse script arguments."""
    parser = ArgumentParser(description="Calculate CPM per contig",
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i", "--samtools_idxstats", nargs='+',  required = True, 
    help = "samtools idxstats file containing contig id, \
    sequence length, number of mapped reads or fragments, \
    number of unmapped reads or fragments")

    parser.add_argument('-f', '--fasta_file',  required = True, 
                        help = 'fasta file containing sequences of contigs.')
    parser.add_argument("-c", "--cutoff_cpm", required = True,
                       help = "Minimum number of reads in a contig")
    parser.add_argument("-s", "--select", 
                       help = "Name of outpout .fa file containing contigs which passed cpm cutoff")
    parser.add_argument("-d", "--discard", 
                       help = "Name of outpout .fa file containing contigs which don't passed cpm cutoff")


    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")

    args = parser.parse_args()
    return args

def combine_idxstat_files(idxstat_dfs):
    """
    Combine multiple idxstat df that have the same contigs.

    Sum the #_mapped_read_segments column over multiple idxstat files that have the same reference sequences.
    """
    
    common_contigs = set(idxstat_dfs[0].index)
    for df in idxstat_dfs[1:]:
        common_contigs = common_contigs.union(set(df.index))
        
    # Convert set to list for pandas indexing
    common_contigs_list = sorted(list(common_contigs))
    
    # Filter dataframes to keep only common contigs
    filtered_dfs = [df.loc[common_contigs_list] for df in idxstat_dfs]
    
    # Combine dataframes
    combined_df = filtered_dfs[0].copy()
    
    for df in filtered_dfs[1:]:
        combined_df['#_mapped_read_segments'] += df['#_mapped_read_segments']
        combined_df['cpm_count'] += df['cpm_count']
        
    return combined_df      
   
def combine_kept_contigs(kept_contigs_list):
    """Combine multiple lists of kept contigs."""
    # Get contigs that appear in all lists (union)
    common_contigs = set(kept_contigs_list[0])
    for contigs in kept_contigs_list[1:]:
        common_contigs = common_contigs.union(set(contigs))
    return pd.Index(list(common_contigs))

def filter_cpm(idxstat_file, cpm_cutoff):
    """
    Calculate and filter byCPM for each contig.
    """
    columns_names = ['reference_sequence_name', 
                    'sequence_length',
                    '#_mapped_read_segments',    
                    '#_unmapped_read-segments']

    idxstat_df = pd.read_csv(idxstat_file, sep ='\t',
                            names = columns_names, 
                            usecols = ['reference_sequence_name', 
                                        'sequence_length',
                                        '#_mapped_read_segments',],
                            comment="*").set_index('reference_sequence_name')
    
    
    logging.info(f'idxstat_file: {idxstat_df}')
    
    
    sum_reads = idxstat_df['#_mapped_read_segments'].sum()
    idxstat_df['cpm_count'] = 1e6 * (idxstat_df['#_mapped_read_segments'] / sum_reads)
    
    # Apply CPM cutoff
    kept_contigs = idxstat_df.loc[idxstat_df["cpm_count"] >= cpm_cutoff].index
    
    logging.info(f'length of idxstat_file: {len(idxstat_df)}')
    logging.info(f'length of kept_contigs: {len(kept_contigs)}')

    return kept_contigs, idxstat_df


def main():
    args = parse_arguments()

    if args.verbose:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG)
        logging.info('Mode verbose ON')

    else:
        logging.basicConfig(format="%(levelname)s: %(message)s")

    cpm_cutoff = float(args.cutoff_cpm)
    
    # logging.info(f'idxstat file: {idxstat_file}')
    idxstat_results = [filter_cpm(idxstat_file, cpm_cutoff) for idxstat_file in args.samtools_idxstats]
    
    # Separate idxstat and kept_contigs
    kept_contigs, idxstat_dfs = zip(*idxstat_results)
        
    logging.info(f'number of kept_contigs: {len(kept_contigs)}')
    logging.info(f'number of idxstat_file: {len(idxstat_dfs)}')
    
    # Combine idxstat dataframes
    combined_idxstat_df = combine_idxstat_files(idxstat_dfs)   
    logging.info(f'length of list_idxstat_df: {len(combined_idxstat_df)}')
    
    combined_kept_contigs = combine_kept_contigs(kept_contigs)
    combined_kept_contigs = set(combined_kept_contigs.tolist())
    logging.info(f'length of combined_kept_contigs: {len(combined_kept_contigs)}')
    
    logging.info(f'{len(combined_kept_contigs)}/{len(combined_idxstat_df)} contigs are kept with a cpm cutoff of {cpm_cutoff}.')
    
    # Write new fasta files with kept and unkept contigs
    with open(args.select, "w") as out_select_handle, open(args.discard, "w") as out_discard_handle:
        
        for contig, seq in pyfastx.Fasta(args.fasta_file, build_index=False):
            if contig in combined_kept_contigs:
                out_select_handle.write(f'>{contig}\n{seq}\n')
            else:
                out_discard_handle.write(f'>{contig}\n{seq}\n')

if __name__ == "__main__":
    main()