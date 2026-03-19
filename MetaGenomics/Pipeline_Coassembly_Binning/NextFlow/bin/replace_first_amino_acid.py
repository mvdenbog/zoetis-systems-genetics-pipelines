#!/usr/bin/env python

"""----------------------------------------------------------------------------
  Script Name: translate_fasta.py
  Description: Corrects the first amino acid of protein sequences in a FASTA 
               file by re-translating the first codon from the corresponding 
               nucleotide sequence. For each protein record, looks up the 
               matching nucleotide record by ID, extracts the first three 
               nucleotides, translates them using the standard genetic code, 
               and replaces the first amino acid of the protein sequence with 
               this newly translated residue. This ensures consistency between 
               nucleotide and protein sequences when the initial translation 
               may have used alternative start codon handling or frame offsets. 
               Outputs a corrected protein FASTA file with all other amino 
               acids preserved unchanged. Requires that every protein ID in 
               the input protein file has a matching nucleotide record.
  Input files: 
    - Nucleotide FASTA (ffn): coding sequences with IDs matching protein file
    - Protein FASTA (faa): translated sequences requiring first-AA correction
  Output files:
    - Corrected protein FASTA (faa): same IDs and descriptions, with first 
      amino acid re-translated from nucleotide first codon
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
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

### Functions
def correct_first_codon(sequence):
    """
    Correct the first codon of the nucleotide sequence.
    
    Args:
    sequence (Seq): Nucleotide sequence.
    
    Returns:
    Seq: Protein sequence for the first codon.
    """
    first_codon = sequence[:3]
    return first_codon.translate()

def replace_first_amino_acid(nucleotide_fasta, protein_fasta, output_fasta):
    nucleotide_records = SeqIO.to_dict(SeqIO.parse(nucleotide_fasta, "fasta"))
    modified_protein_records = []

    # Read the protein sequences from the FASTA file
    for protein_record in SeqIO.parse(protein_fasta, "fasta"):
        if protein_record.id in nucleotide_records:
            nucleotide_seq = nucleotide_records[protein_record.id].seq
            # Translate the first codon into amino acid
            first_amino_acid = correct_first_codon(nucleotide_seq)
            # Replace the first amino acid of the protein sequence
            modified_protein_seq = first_amino_acid + protein_record.seq[1:]
            
            # Create a new record for the modified protein
            modified_protein_record = SeqRecord(
                Seq(modified_protein_seq),
                id=protein_record.id,
                description=protein_record.description
            )
            
            # Add the modified protein record to the list
            modified_protein_records.append(modified_protein_record)
    
    # Write the modified protein sequences to the output FASTA file
    SeqIO.write(modified_protein_records, output_fasta, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: replace_first_amino_acid.py <nucleotide_fasta> <protein_fasta> <output_fasta>")
        sys.exit(1)

    nucleotide_fasta = sys.argv[1]
    protein_fasta = sys.argv[2]
    output_fasta = sys.argv[3]

    replace_first_amino_acid(nucleotide_fasta, protein_fasta, output_fasta)
    print(f"Modified protein sequences written in {output_fasta}")