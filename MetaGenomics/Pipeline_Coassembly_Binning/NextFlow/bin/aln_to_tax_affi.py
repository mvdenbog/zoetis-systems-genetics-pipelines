#!/usr/bin/env python3

"""----------------------------------------------------------------------------
  Script Name: aln_to_tax_affi.py
  Description: Assigns taxonomic labels to metagenomic contigs and their 
               constituent proteins based on Diamond/BLAST alignment results 
               (m8 format). Uses a weighted scoring system that normalizes 
               alignment scores by query coverage and identity, then propagates 
               taxonomic weights up the NCBI taxonomy tree. Consensus taxonomy 
               is determined per protein and per contig using a minimum fraction 
               threshold (default: 90%) of cumulative weight at each rank 
               (kingdom to species). Supports filtering by minimum identity, 
               coverage, and alignment score range. Optionally outputs top 
               candidate taxa per rank with weights for interpretability.
  Input files: 
    - Diamond/BLAST alignment file (.m8 format with qlen column)
    - Accession-to-taxid mapping file (gzipped or plain TSV)
    - Extracted NCBI taxdump directory (nodes.dmp, names.dmp, merged.dmp, taxidlineage.dmp)
    - Optional: query length file for plotting statistics
  Output files:
    - <output>.pergene.tsv: per-protein taxonomic assignments
    - <output>.percontig.tsv: per-contig consensus taxonomic assignments  
    - <output>.warn.tsv: warnings for accessions missing taxid mapping or taxonomy
    - Optional: top_taxons_per_contig.tsv / _verbose.tsv with weighted candidates
    - Optional: PDF plots of assignment statistics (if query_length_file provided)
  Created By:  Mathias Vandenbogaert
  Date:        2026-05-06
-------------------------------------------------------------------------------
"""

# Metadata
__author__ = 'Mathias Vandenbogaert'
__version__ = '0.1'
__email__ = 'mathias.vandenbogaert@lizard.bio'
__status__ = 'dev'


# Modules importation

import gzip
import re
import os
import operator
from collections import defaultdict
from collections import Counter
import csv
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
from matplotlib import pyplot
from collections import OrderedDict
import pandas as pd

# Variables

# These are identities normalized with query coverage:
RANKS_TO_MIN_SCORE = {'kingdom': 0.4,
                      'phylum': 0.5,
                      'class': 0.6,
                      'order': 0.7,
                      'family': 0.8,
                      'genus': 0.9,
                      'species': 0.95}

# Fraction of weights needed to assign at a specific level,
# a measure of concensus at that level.
MIN_FRACTION = 0.9


def parse_arguments():
    """Parse script arguments."""
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument("-b", "--aln_input_file", required=True,
                        help="file with blast/diamond matches expected format m8 \
    \nqueryId, subjectId, percIdentity, alnLength, mismatchCount, gapOpenCount,\
    queryStart, queryEnd, subjectStart, subjectEnd, eVal, bitScore")

    parser.add_argument('-a', '--acc_taxaid_mapping_file', required=True,
                        help="mapping from accession to taxaid gzipped")

    parser.add_argument('-t', '--taxonomy', required=True,
                        help="path of taxdump.tar.gz extracted directory")

    parser.add_argument('-o', '--output_file', type=str,
                        default="taxonomyassignation", help=("string specifying output file"))

    parser.add_argument('-i', '--min_identity', default=60,
                        help="percentage of identity")

    parser.add_argument('-c', '--min_coverage', default=70,
                        help="percentage of coverage")

    parser.add_argument('--top', default=10, type=int,
                        help="Keep diamond alignments within this percentage range of top alignment score")

    parser.add_argument('--keep_only_best_aln',
                        help="Keep only diamond alignments with top alignment score. (overrides --top)", action="store_true")

    parser.add_argument('--write_top_taxons',
                    help="""Write top taxons per contig for each rank
                    with their weigth associated in 'top_taxons_per_contig.tsv'.
                    Can be helpful to understand the affiliations made.""",
                    action="store_true")

    parser.add_argument('--write_top_taxons_verbose',
                    help="""Write top taxons per contig for each rank
                    with their weigth associated in a verbose mode 'top_taxons_per_contig_verbose.tsv'.
                    Can be helpful to understand the affiliations made.""",
                    action="store_true")

    parser.add_argument("--query_length_file",
                        help="tab delimited file of query lengths")

    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")

    args = parser.parse_args()

    return args


################################################
# Functions for taxonomy
################################################

def load_taxonomy(taxdump_dir, main_ranks, taxids_selection):

    logging.info(f'Load taxonomy information for {len(taxids_selection)} taxids.')

    nodes_file = os.path.join(taxdump_dir, "nodes.dmp")
    taxidlineage_file = os.path.join(taxdump_dir, "taxidlineage.dmp")
    merged_file = os.path.join(taxdump_dir, "merged.dmp")

    names_file = os.path.join(taxdump_dir, "names.dmp")

    merged_taxid = replace_merged_taxid(taxids_selection, merged_file)
    logging.info(f'{len(merged_taxid)} taxids have been merged into another taxids: {merged_taxid}')

    taxid2lineage_whole_db = get_all_taxid_lineage(taxidlineage_file)

    all_taxids = {taxid for leaf_taxid, lineage in taxid2lineage_whole_db.items()
                  for taxid in lineage if leaf_taxid in taxids_selection}

    all_taxids.add(1)
    taxid2lineage = {taxid: lineage for taxid,
                     lineage in taxid2lineage_whole_db.items() if taxid in all_taxids}

    taxid2rank = get_taxid_rank(all_taxids, nodes_file)

    taxid2rankedlineage = {taxid: get_ranked_lineage(
        lineage, taxid2rank, main_ranks) for taxid, lineage in taxid2lineage.items()}

    taxid2name = get_taxid2name(all_taxids, names_file)

    logging.info(f'Load taxid information done. {len(taxids_selection - set(taxid2rankedlineage))} taxid has not been found in taxdump files')
    return taxid2rankedlineage, taxid2name, taxid2rank, merged_taxid


def get_taxid2name(taxids, names_file):
    taxid2name = {"None": "None"}
    with open(names_file, "r") as name_file:
        for line in name_file:
            line = line.rstrip().replace("\t", "")
            tab = line.split("|")
            if int(tab[0]) in taxids and tab[3] == "scientific name":
                tax_id, name = int(tab[0]), tab[1]
                taxid2name[tax_id] = name
    return taxid2name


def replace_merged_taxid(taxids, merged_file):
    merged_taxid = set()
    with open(merged_file) as fl:
        for i, l in enumerate(fl):
            old_taxid, new_taxid = l.rstrip().replace('\t|', '').split('\t')

            if int(old_taxid) in taxids:
                taxids.remove(int(old_taxid))
                taxids.add(int(new_taxid))

                merged_taxid.add(int(old_taxid))
    return merged_taxid


def get_ranked_lineage(taxids, taxid2rank, ranks_to_keep):
    rank2taxid = {taxid2rank[taxid]: taxid for taxid in taxids}
    ranked_lineage = []
    for rank_to_keep in ranks_to_keep:
        try:
            ranked_lineage.append(rank2taxid[rank_to_keep])
        except KeyError:
            ranked_lineage.append('None')
    return ranked_lineage


def get_taxid_rank(taxids, nodes):
    taxid2rank = {}
    with open(nodes) as fl:
        for i, l in enumerate(fl):
            node_infos = l.rstrip().replace('\t|', ' ').split('\t')
            taxid = node_infos[0].strip()
            if int(taxid) in taxids:
                rank = node_infos[2].strip()
                taxid2rank[int(taxid)] = rank

    return taxid2rank


def get_all_taxid_lineage(taxid_lineage):

    taxid2lineage = {}
    with open(taxid_lineage) as fl:
        for i, l in enumerate(fl):
            taxid, taxid_lineage_str = l.rstrip().replace('\t|', ' ').split('\t')
            taxid_lineage = [int(taxid) for taxid in taxid_lineage_str.strip().split(' ') if taxid]
            taxid_lineage.append(int(taxid))
            taxid2lineage[int(taxid)] = taxid_lineage

    return taxid2lineage

################################################
# END Functions for taxonomy
################################################


def read_blast_input(blastinputfile, min_identity, min_coverage, top_aln):

    logging.info(f'Parsing blast result file {blastinputfile}...')

    matches = defaultdict(list)
    # accs = Counter()
    min_score_per_query = defaultdict(int)

    with open(blastinputfile) as blast_handler:

        reader = csv.DictReader(blast_handler, delimiter='\t')

        for aln in reader:

            if aln['sseqid'].startswith("gi|"):
                m = re.search(r"gi\|.*?\|.*\|(.*)\|", aln['sseqid'])
                acc = m.group(1)
            else:
                acc = aln['sseqid']

            qlen = int(aln['qlen'])
            qseqid = aln['qseqid']

            aln_qlen = abs(int(aln['qend']) - int(aln['qstart'])) + 1
            qcov = 100 * float(aln_qlen) / qlen
            pident = float(aln['pident'])

            if pident >= min_identity and qcov >= min_coverage:

                score = qcov/100 * pident/100

                assert score <= 1.0

                # when top = 0, the  min_score_per_query is the best score
                min_score_per_query[qseqid] = max(
                    min_score_per_query[qseqid], score*(1 - top_aln/100))

                if min_score_per_query[qseqid] <= score:
                    matches[qseqid].append({'sseqid': acc, 'score': score})

    logging.info(f'Parsing blast result file {blastinputfile} is completed. {len(matches)} proteins found with diamond hits.')

    matches_filtered = {}
    for qseqid,  infos in matches.items():
        matches_filtered[qseqid] = [
            info for info in infos if info['score'] >= min_score_per_query[qseqid]]

    return matches_filtered


def parse_accession2taxid(accs, mapping_file):

    logging.info(f'Parsing accession2taxid file {mapping_file} to retrieve taxid of {len(accs)} protein accessions.')

    counter = 0
    total_prot = len(accs)
    proper_open = gzip.open if mapping_file.endswith('.gz') else open
    mappings = {acc: None for acc in accs}
    with proper_open(mapping_file, 'rt') as mapping_fh:
        for line in mapping_fh:
            acc_ver, taxid = line.split("\t")
            # Only add taxids for the given acc
            if acc_ver in mappings:
                counter += 1
                if counter % (total_prot/10) == 0:
                    logging.info(f'{100 * counter/total_prot:.0f}% DONE')
                mappings[acc_ver] = int(taxid)

    logging.info(f'Parsing accession2taxid is completed. {len([t for t in mappings.values() if not t])}/{len(accs)} accessions have no taxid associated.')

    return mappings


def group_by_contig(matches, contig_pattern):
    contig2matches = defaultdict(dict)

    for prot_id, hits in matches.items():
        contig = contig_pattern.match(prot_id).group(1)
        contig2matches[contig][prot_id] = hits
    return contig2matches


def collate_protein_hits(sorted_hits, main_ranks, taxid2lineage, accession2taxid):
    """Collatte protein hits."""

    taxids_already_processed = set()

    collate_hits = {rank: Counter() for rank in main_ranks}

    # For each hit, retrieve taxon id and compute weight in lineage
    for hit in sorted_hits:
        logging.debug('====HIT====')
        protein_hit = hit['sseqid']
        score = hit['score']
        logging.debug(f'{protein_hit}, {score}')

        hit_taxid = accession2taxid[protein_hit]

        if hit_taxid is None:
            logging.debug(f'    {protein_hit} accession has no corresponding taxid in accession2taxid file')
            #  protein hit has not been found in accession2taxid
            continue

        if hit_taxid in taxids_already_processed:
            logging.debug(f'    {hit_taxid} already processed')
            # Only add the best hit per species
            continue

        taxids_already_processed.add(hit_taxid)

        if hit_taxid in taxid2lineage:
            logging.debug('    Hit has a taxo')
            hit_taxonomy = taxid2lineage[hit_taxid]

            logging.debug(f'    Protein taxo: {hit_taxonomy}')
            for rank, rank_taxid in zip(main_ranks, hit_taxonomy):
                if rank_taxid == "None":
                    continue

                weight = (score - RANKS_TO_MIN_SCORE[rank]) / (1.0 - RANKS_TO_MIN_SCORE[rank])
                weight = max(weight, 0.0)

                logging.debug(f"    R={rank}, T={rank_taxid}, W={weight}")

                # could put a transform in here
                if weight > 0:
                    collate_hits[rank][rank_taxid] += weight

        else:
            logging.debug(f'        Hit taxid {hit_taxid} is not found in taxo')
    return collate_hits


def get_taxid_consensus(collate_table, main_ranks):

    for rank in main_ranks[::-1]:
        collate = collate_table[rank]
        if not collate:
            continue
        dWeight = sum(collate.values())
        sortCollate = sorted(list(collate.items()), key=operator.itemgetter(1), reverse=True)
        logging.debug(f"{rank}, {sortCollate}, sum score {dWeight}")
        best_taxid, best_taxid_score = sortCollate[0]

        if len(collate) > 0 and dWeight > 0.0:
            dP = float(best_taxid_score) / dWeight
            if dP > MIN_FRACTION:
                logging.debug(f'-->dP OK {best_taxid}')
                return best_taxid
                # (fullnamelineage_text, fullnamelineage_ids) = d_taxonomy[str(sortCollate[0][0])].lineage_main_level()
                # tax_id_keep = str(sortCollate[0][0])
                # return (tax_id_keep, fullnamelineage_text, fullnamelineage_ids)
    return 1  # (1," Unable to find taxonomy consensus",1)


def get_top_taxid(collate_table, main_ranks):
    top_taxons_per_rank = {}

    for rank in main_ranks[::-1]:
        collate = collate_table[rank]
        if not collate:
            continue

        dWeight = sum(collate.values())
        sortCollate = sorted(list(collate.items()), key=operator.itemgetter(1), reverse=True)
        logging.debug(f"{rank}, {sortCollate}, sum score {dWeight}")

        top_taxids = [(taxid, taxid_score/dWeight) for taxid, taxid_score in sortCollate if taxid_score/dWeight > 0.01]
        top_taxons_per_rank[rank] = top_taxids

    return top_taxons_per_rank  # (1," Unable to find taxonomy consensus",1)


def add_collate_hits(main_collate_hits, collate_hits_to_add):
    for rank in main_collate_hits:
        logging.debug('RANK {rank}')
        for taxid in collate_hits_to_add[rank]:
            logging.debug(f'  {taxid} {collate_hits_to_add[rank][taxid]}')
            main_collate_hits[rank][taxid] += collate_hits_to_add[rank][taxid]


def get_affilaition_line(contig, taxid, taxid2rankedlineage, taxid2name):

    if taxid == 1:
        return f'{contig}\t{taxid}\t Unable to find taxonomy consensus\t{taxid}\n'

    taxid_lineage = taxid2rankedlineage[taxid]
    named_lineage_text = '; '.join([taxid2name[taxid] for taxid in taxid_lineage])

    taxid_lineage_str = '; '.join([str(t) for t in taxid_lineage])

    contig_affi_line = f'{contig}\t{taxid}\t{named_lineage_text}\t{taxid_lineage_str}\n'
    return contig_affi_line


def plot_taxonomic_assignment(output_name, count_genealogy,  count_genealogy_contig, nb_total_prot, nb_prot_annotated, nb_prot_assigned):
    # graphs
    try:
        os.makedirs("graphs")
    except OSError:
        if not os.path.isdir("graphs"):
            Raise

    # Sort dictionaries
    count_genealogy_ord = OrderedDict(sorted(count_genealogy.items(), key=lambda t: t[0]))
    count_genealogy_contig_ord = OrderedDict(
        sorted(count_genealogy_contig.items(), key=lambda t: t[0]))
    # Figures
    pyplot.bar(range(len(count_genealogy_ord.values())), count_genealogy_ord.values())
    pyplot.xticks(range(len(count_genealogy_ord.values())), count_genealogy_ord.keys())
    pyplot.xlabel("Taxonomy level")
    pyplot.ylabel("Number of proteins")
    pyplot.title(output_name + " number of proteins at different taxonomy levels")
    pyplot.savefig("graphs/" + output_name + "_prot_taxonomy_level.pdf")
    pyplot.close()
    pyplot.bar(range(len(count_genealogy_contig_ord.values())), count_genealogy_contig_ord.values())
    pyplot.xticks(range(len(count_genealogy_contig_ord.values())),
                  count_genealogy_contig_ord.keys())
    pyplot.xlabel("Taxonomy level")
    pyplot.ylabel("Number of contigs")
    pyplot.title(output_name + " number of contigs at different taxonomy levels")
    pyplot.savefig("graphs/" + output_name + "_contig_taxonomy_level.pdf")
    pyplot.close()

    list_graphs = [nb_total_prot, nb_prot_annotated, nb_prot_assigned]
    pyplot.bar(range(len(list_graphs)), list_graphs)
    pyplot.xticks(range(len(list_graphs)), ["Total", "Annotated", "Assigned"])
    pyplot.ylabel("Number of proteins")
    pyplot.title(output_name + " number of annotated and assigned proteins")
    pyplot.savefig("graphs/" + output_name + "_nb_prot_annotated_and_assigned.pdf")
    pyplot.close()


def get_top_taxons_info(contig, top_taxons_per_rank, taxid2name):

    info = {'contig': contig, }
    for rank, top_taxids in top_taxons_per_rank.items():
        top_affi_by_rank = []
        for taxid, score in top_taxids:
            taxname_and_weigth = f"{taxid2name[taxid]} ({100*score:.1f})"
            top_affi_by_rank.append(taxname_and_weigth)
        info[rank] = ';'.join(top_affi_by_rank)
    return info


def get_top_taxons_info_verbose(contig, top_taxons_per_rank, taxid2name, taxid2rankedlineage):
    list_info = []
    for rank, top_taxids in top_taxons_per_rank.items():

        for taxid, score in top_taxids:

            lienage = [taxid2name[tid] for tid in taxid2rankedlineage[taxid] if tid!="None"]

            info= {"contig": contig,
                "rank": rank,
                "lineage": ';'.join(lienage),
                'taxon': taxid2name[taxid],
                "score": round(100 * score, 1)
            }
            list_info.append(info)
    return list_info


def main():

    args = parse_arguments()

    if args.verbose:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
        logging.info('Mode verbose ON')

    else:
        logging.basicConfig(format="%(levelname)s: %(message)s")

    aln_input_file = args.aln_input_file

    min_identity = args.min_identity
    min_coverage = args.min_coverage
    top_aln = args.top
    keep_only_best_aln_flag = args.keep_only_best_aln
    query_length_file = args.query_length_file

    top_taxon_outfile = 'top_taxons_per_contig.tsv'
    top_taxon_verbose_outfile = 'top_taxons_per_contig_verbose.tsv'

    if keep_only_best_aln_flag:
        top_aln = 0

    acc_taxaid_mapping_file = args.acc_taxaid_mapping_file
    taxdump_dir = args.taxonomy

    output_name = args.output_file

    prot_prefix = 'CDS_'

    main_ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

    matches = read_blast_input(aln_input_file, min_identity, min_coverage, top_aln)

    accessions = {hit['sseqid'] for hits in matches.values() for hit in hits}
    accession2taxid = parse_accession2taxid(accessions, acc_taxaid_mapping_file)

    accession_taxids = {t for t in accession2taxid.values() if t}
    taxid2rankedlineage, taxid2name, taxid2rank, merged_taxid = load_taxonomy(
        taxdump_dir, main_ranks, accession_taxids)

    accessions_unfound_in_mapping = {acc for acc, taxid in accession2taxid.items() if taxid is None}

    taxids_found_in_taxonomy = (set(taxid2rankedlineage) | merged_taxid)
    acessions_with_unfound_taxid = {(acc) for acc, taxid in accession2taxid.items(
    ) if taxid not in taxids_found_in_taxonomy and taxid is not None}

    logging.info(f'{len(accessions_unfound_in_mapping)} protein accessions have not been found in accession2taxid file')
    logging.info(f'{len(acessions_with_unfound_taxid)} protein accessions found in accession2taxid have their taxid not found in taxdump')

    re_contig = re.compile('(.*)\.' + prot_prefix)

    count_rank_affiliation_protein = Counter()
    count_rank_affiliation_contig = Counter()
    nb_prot_assigned = 0

    contig2matches = group_by_contig(matches, re_contig)
    top_taxon_infos = []
    top_taxon_infos_verbose = []

    with open(output_name + ".pergene.tsv", "w") as out_protein, \
            open(output_name + ".percontig.tsv", "w") as out_contig, \
            open(output_name + ".warn.tsv", "w") as outdisc:

        # Write header
        out_protein.write("#prot_id\tconsensus_tax_id\tconsensus_lineage\ttax_id_by_level\n")
        out_contig.write("#contig\tconsensus_tax_id\tconsensus_lineage\ttax_id_by_level\n")
        outdisc.write("#prot_id\tlist nr hit not found in taxo\n")

        for contig, contig_matches in contig2matches.items():
            logging.info(contig)
            contig_collate_hits = {rank: Counter() for rank in main_ranks}

            for protein_id, hits in contig_matches.items():

                sorted_hits = sorted(hits, key=lambda x: x["score"], reverse=True)

                # manage protein affiliation
                protein_collate_hits = collate_protein_hits(
                    hits, main_ranks, taxid2rankedlineage, accession2taxid)
                consensual_protein_taxid = get_taxid_consensus(protein_collate_hits, main_ranks)

                count_rank_affiliation_protein[taxid2rank[consensual_protein_taxid]] += 1

                protein_affi_line = get_affilaition_line(
                    protein_id, consensual_protein_taxid, taxid2rankedlineage, taxid2name)

                out_protein.write(protein_affi_line)
                nb_prot_assigned += 1

                logging.debug(f'PROTEIN {consensual_protein_taxid}, {taxid2rankedlineage[consensual_protein_taxid]}')

                # manage contig affiliation
                protein_collate_hits = collate_protein_hits(
                    sorted_hits, main_ranks, taxid2rankedlineage, accession2taxid)

                add_collate_hits(contig_collate_hits, protein_collate_hits)

                hit_accessions = {hit['sseqid'] for hit in hits if hit['sseqid']}
                acessions_with_unfound_taxid_prot = acessions_with_unfound_taxid & hit_accessions
                accessions_unfound_in_mapping_prot = accessions_unfound_in_mapping & hit_accessions

                if acessions_with_unfound_taxid_prot:
                    outdisc.write(f"{protein_id}\tNo taxid in taxdump\t{','.join(sorted(acessions_with_unfound_taxid_prot))}\n")

                if accessions_unfound_in_mapping_prot:
                    outdisc.write(f"{protein_id}\tNo protid correspondance file\t{','.join(sorted(accessions_unfound_in_mapping_prot))}\n")

            consensual_contig_taxid = get_taxid_consensus(contig_collate_hits, main_ranks)
            if args.write_top_taxons or args.write_top_taxons_verbose:
                top_taxons_per_rank = get_top_taxid(contig_collate_hits, main_ranks)

            if args.write_top_taxons:
                top_taxon_info = get_top_taxons_info(contig, top_taxons_per_rank, taxid2name)
                top_taxon_infos.append(top_taxon_info)

            if args.write_top_taxons_verbose:
                top_taxon_infos_verbose += get_top_taxons_info_verbose(contig, top_taxons_per_rank, taxid2name, taxid2rankedlineage)

            count_rank_affiliation_contig[taxid2rank[consensual_contig_taxid]] += 1
            logging.debug(f'CONTIG {consensual_contig_taxid}, {taxid2rankedlineage[consensual_contig_taxid]}')
            contig_affi_line = get_affilaition_line(
                contig, consensual_contig_taxid, taxid2rankedlineage, taxid2name)
            out_contig.write(contig_affi_line)
            logging.debug(contig_affi_line)

    if query_length_file:
        logging.debug("Plot taxonomic affiliation using protein lengths.")
        with open(query_length_file) as fl:
            nb_total_prot = len([line for line in fl])

        nb_prot_annotated = len(matches)
        plot_taxonomic_assignment(
            output_name, count_rank_affiliation_protein,  count_rank_affiliation_contig, nb_total_prot, nb_prot_annotated, nb_prot_assigned)

    if args.write_top_taxons:
        top_taxon_columns = ['contig', ] + main_ranks
        pd.DataFrame(top_taxon_infos, columns=top_taxon_columns).to_csv(top_taxon_outfile, sep='\t', index=False)

    if args.write_top_taxons_verbose:
        pd.DataFrame(top_taxon_infos_verbose).to_csv(top_taxon_verbose_outfile, sep='\t', index=False)


if __name__ == '__main__':
    main()