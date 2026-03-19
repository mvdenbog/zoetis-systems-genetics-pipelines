#!/usr/bin/env python

"""----------------------------------------------------------------------------
  Script Name: db_info_summary.py
  Description: Summarizes metadata for a list of biological databases specified 
               in an input text file. For each database entry (name and path), 
               it retrieves the disk size (using 'du -sh'), last modification 
               date (using 'stat'), and resolved absolute path. Special handling 
               is applied: if the path points to a 'nodes.dmp' file, the size 
               of the parent directory is reported; if the database name is 
               'Host_genome', the total sequence count (lines starting with '>') 
               is appended to the size string. Supports gzipped genome files. 
               Outputs a tab-separated summary to stdout including the original 
               input line for reference.
  Input files: 
    - Text file listing databases (format: 'name path' per line)
    - Referenced database files/directories (must exist on filesystem)
  Output: 
    - Stdout: Header line followed by original input lines and corresponding 
      metadata (Name, Size, Last modification, Path) tab-separated.
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
import subprocess
import os.path
import gzip


def get_genome_seq_count(genome_path):
    """
    Count the number of sequences in a genome file (fasta format).
    Supports gzipped files (.gz).
    """
    seq_count = 0

    proper_open = gzip.open if genome_path.endswith('.gz') else open
    with proper_open(genome_path, "rt") as fh:
        for l in fh:
            if l.startswith('>'):
                seq_count += 1
    return seq_count


def info_db(db_info_file):
    """
    Extract metadata for a single database entry.
    
    Args:
        db_info_file: String containing 'name path' (space-separated).
        
    Returns:
        Tuple: (name, size, modif_date, path)
    """
    name = db_info_file.split()[0]
    path = os.path.realpath(db_info_file.split()[1])

    if db_info_file.split()[1] == 'nodes.dmp':
        path = os.path.dirname(path)

    # get db size 
    process = subprocess.run(['du', '-sh', path], stdout=subprocess.PIPE)
    size = process.stdout.split()[0].decode('utf-8')

    if (name == "Host_genome"):
        size = f"{size} ({get_genome_seq_count(path)} seq)"

    # get date of last modifaction 
    process = subprocess.run(['stat', '-c %y', path], stdout=subprocess.PIPE)
    modif_date = process.stdout.split()[0].decode('utf-8')

    return name, size, modif_date, path


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--all_db', required=True, help=
    'text file with details of all databank used ')
    args = parser.parse_args()

    print("DB (folder or file)\tSize\tLast modification\tPath")

    with open(args.all_db, "r") as list_db:
        lines = list_db.readlines()
        for line in lines:
            print(line)
            list_info = info_db(line.strip())
            print(*list_info, sep="\t")


if __name__ == "__main__":
    main()