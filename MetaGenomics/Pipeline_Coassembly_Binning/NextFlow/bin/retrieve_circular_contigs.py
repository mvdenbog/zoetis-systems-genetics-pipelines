#!/usr/bin/env python3

"""----------------------------------------------------------------------------
  Script Name: extract_circular_contigs.py
  Description: Extracts circular contigs from long-read metagenomic assemblies 
               produced by different assemblers (metaFlye, hifiasm-meta, 
               metaMDBG). Circular contig identification is assembler-specific: 
               (1) metaFlye: parses assembly_info.txt file and selects contigs 
               with 'Y' or '+' in the 'circ.' column; (2) hifiasm-meta: 
               identifies contigs with 'c' suffix in the FASTA description; 
               (3) metaMDBG: same 'c' suffix detection but handles gzipped 
               FASTA input. Each circular contig is written to a separate 
               FASTA file (bin_<N>.fa) in the output directory, enabling 
               downstream plasmid analysis, circular genome validation, or 
               quality assessment of assembly completeness. Supports both 
               plain and gzipped FASTA inputs.
  Input files: 
    - Assembly FASTA file: contig sequences from metaFlye, hifiasm-meta, or 
      metaMDBG (gzipped for metaMDBG)
    - For metaFlye only: assembly_info.txt file with circ. column indicating 
      circularity status
  Output files:
    - Output directory (default: circular_contigs/) containing one FASTA file 
      per circular contig (bin_0.fa, bin_1.fa, etc.)
  Created By:  Mathias Vandenbogaert
  Date:        2026-05-06
-------------------------------------------------------------------------------
"""

# Metadata
__author__ = 'Mathias Vandenbogaert'
__version__ = '0.1'
__email__ = 'mathias.vandenbogaert@lizard.bio'
__status__ = 'dev'

## Imports
from Bio import SeqIO
import pandas as pd
import argparse
import gzip
import os

## Functions

def parse_metaflye(input_fasta, info_file, outdir):
    """
    metaflye fasta contains no information.
    This is contained in the assembly_info.txt file.
    """
    df = pd.read_csv(info_file, sep='\t')
    # this file may have changed symbols for yes, possibilities include a Y or +
    circular_contigs = list(df[df['circ.'].isin(['Y', '+'])]["#seq_name"])

    for count, rec in enumerate(SeqIO.parse(input_fasta, "fasta")):
        if rec.id in circular_contigs:
            outfile = os.path.join(outdir, f"bin_{count}.fa")
            with open(outfile, 'wt') as outfl:
                outfl.write(rec.format("fasta"))

    # with open(circular_fasta, 'wt') as fh_out:
    #     for rec in SeqIO.parse(input_fasta, "fasta"):
    #         if rec.id in circular_contigs:
    #             fh_out.write(rec.format("fasta"))

def parse_hifiasm(input_fasta, outdir):
    """
    hifiasm-meta circular contig names will have an "c" suffix.
    """
    for count, rec in enumerate(SeqIO.parse(input_fasta, "fasta")):
        if rec.description.endswith('c'):
            outfile = os.path.join(outdir, f"bin_{count}.fa")
            with open(outfile, 'wt') as outfl:
                outfl.write(rec.format("fasta"))

    # with open(circular_fasta, 'wt') as fh_out:
    #     for rec in SeqIO.parse(input_fasta, "fasta"):
    #             if rec.id.endswith('c'):
    #                 fh_out.write(rec.format("fasta"))

def parse_metamdbg(input_fasta, outdir):
    """
    metaMDBG circular contig names will have an "c" suffix but they are in
    fasta.gz.
    """
    with gzip.open(input_fasta, "rt") as inp_f:
        for count, rec in enumerate(SeqIO.parse(inp_f, "fasta")):
            if rec.description.endswith('c'):
                outfile = os.path.join(outdir, f"bin_{count}.fa")
                with open(outfile, 'wt') as outfl:
                    outfl.write(rec.format("fasta"))

    # with open(circular_fasta, 'wt') as fh_out:
    #     for rec in SeqIO.parse(input_fasta, "fasta"):
    #             if rec.id.endswith('c'):
    #  

def main():
    # Manage parameters
    parser = argparse.ArgumentParser(description='Retrieve circular contigs from long-read assemblies.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-a', '--assembler', required=True, choices=['metaflye', 'hifiasm-meta', 'metamdbg'], help='Assembler where the assemblies come from.')
    group_input.add_argument('-f', '--input-fasta', required=True, help='The path of assembly fasta file.')
    
    group_input_otu_table = parser.add_argument_group(' Metaflye ')
    group_input_otu_table.add_argument('-i','--info-file', default=None, help="The path of metaflye asembly info file.")
    # output
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o','--outdir', default='circular_contigs', help="The path of circular contigs directory. [Default: %(default)s]")

    args = parser.parse_args()

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # Check for inputs
    if args.assembler == "metaflye":
        if args.info_file is None:
            parser.error("\n\n#ERROR : --info-file is required with metaflye assembler.")
        else:
            parse_metaflye(args.input_fasta, args.info_file, outdir)

    if args.assembler == "hifiasm-meta":
        parse_hifiasm(args.input_fasta, outdir)
    
    elif args.assembler == "metamdbg":
        parse_metamdbg(args.input_fasta, outdir)


if __name__ == '__main__':
    main()