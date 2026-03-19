#!/usr/bin/env python3

"""----------------------------------------------------------------------------
  Script Name: filter_diamond_hits.py
  Description: Filters Diamond/BLAST alignment results (m8 format) to retain 
               only the best hit(s) per query sequence based on highest bitscore. 
               Applies quality filters: hits must meet minimum percent identity 
               (default: 60%) and minimum query coverage (default: 70%, calculated 
               as (qend-qstart+1)/qlen * 100). For each query, all hits passing 
               these thresholds are considered, then only those with the maximum 
               bitscore are written to output. Supports ties (multiple hits with 
               identical best bitscore). Requires input file to be sorted by 
               query ID (qseqid) for efficient streaming processing. Outputs a 
               filtered m8-format TSV file with header preserved. Reports summary 
               statistics: number/percentage of queries with no qualifying hits, 
               and total best hits written.
  Input files: 
    - Diamond/BLAST m8 output file (TSV with header: qseqid, sseqid, pident, 
      length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, 
      qlen, slen, stitle); must be sorted by qseqid
  Output files:
    - Filtered m8 TSV file containing only best bitscore hits per query that 
      pass identity and coverage thresholds (default: best_hit.tsv)
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
import csv


def get_hits_with_highest_bitscore(hits):
    highest_bitscore = max([float(hit['bitscore']) for hit in hits])
    return [hit for hit in hits if float(hit['bitscore']) == highest_bitscore]


def get_all_hits_per_query(blast_result_file):
    # Assertion: Hit are already sorted by query in diamond output.
    # Both commands should output the same number of line:
    # cut -f1 blast_result_file | uniq | wc -l
    # cut -f1 blast_result_file | sort | uniq | wc -l

    with open(blast_result_file) as in_fl:

        result_reader = csv.DictReader(in_fl, delimiter='\t')

        query_ids_processed = []

        current_query_id = None
        hits = []

        for hit in result_reader:

            if not current_query_id:
                current_query_id = hit['qseqid']

            if current_query_id and current_query_id != hit['qseqid']:
                yield hits
                hits = []

                current_query_id = hit['qseqid']
                assert current_query_id not in query_ids_processed, f"Queries are not sorted in blast result. Query {current_query_id} is found in different part of the file."

                query_ids_processed.append(current_query_id)

            hits.append(hit)

        if current_query_id:
            yield hits


def is_identity_and_coverage_ok(hit, min_identity, min_coverage):

    qcovhsp = (int(hit["qend"]) - int(hit["qstart"]) + 1) / int(hit['qlen']) * 100
    return float(hit['pident']) >= min_identity and qcovhsp >= min_coverage



def parse_arguments():
    """Parse script arguments."""
    parser = ArgumentParser(description="Filter Diamond hits by identity, coverage, and best bitscore.",
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('aln_input_file',
                        help="File with blast/diamond matches expected format m8 \
    \nqseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore	qlen	slen	stitle")

    parser.add_argument('-o', '--output_file', type=str,
                        default="best_hit.tsv", help=("string specifying output file path"))

    parser.add_argument('-i', '--min_identity', default=60, type=float,
                        help="percentage of identity")

    parser.add_argument('-c', '--min_coverage', default=70, type=float,
                        help="percentage of coverage")

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

    blast_result = args.aln_input_file
    outfile = args.output_file
    min_coverage = args.min_coverage
    min_identity = args.min_identity
    best_hit_count = 0
    query_count_with_low_hit = 0
    with open(blast_result) as f:
        header = f.readline().rstrip().split('\t')
    with open(outfile, 'w') as out_fl:
        
        writer = csv.DictWriter(out_fl, delimiter='\t',fieldnames=header)
        writer.writeheader()
        query_i=0
        for query_i, query_hits in enumerate(get_all_hits_per_query(blast_result)):

            if query_i % 10000 == 0:
                logging.info(f'{query_i} queries processed... ')

            correct_hits = [hit for hit in query_hits if is_identity_and_coverage_ok(
                hit, min_identity, min_coverage)]

            if not correct_hits:
                query_count_with_low_hit += 1
                continue

            best_hits = get_hits_with_highest_bitscore(correct_hits)
            for best_hit in best_hits:
                best_hit_count += 1
                writer.writerow(best_hit)

    logging.info(f'{query_count_with_low_hit} queries ({100*query_count_with_low_hit/(query_i+1):.2f}%) have low hits that do not pass identity ({min_identity}%) or coverage ({min_coverage}%) thresholds')
    logging.info(f'{best_hit_count} best hits of {query_i+1 - query_count_with_low_hit } queries have been written in {outfile}.')


if __name__ == '__main__':
    main()