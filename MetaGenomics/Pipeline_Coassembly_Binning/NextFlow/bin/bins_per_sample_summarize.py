#!/usr/bin/env python

"""----------------------------------------------------------------------------
  Script Name: bins_per_sample_summarize.py
  Description: Generates a comprehensive abundance table of metagenomic bins 
               across multiple samples by aggregating samtools coverage and 
               flagstat outputs. Integrates taxonomic classifications from 
               GTDB-Tk and quality metrics (completeness, contamination, length, 
               N50) from DRep/CheckM2. Outputs: (1) a unified TSV table with 
               per-bin abundances (read counts, mean depth) per sample plus 
               taxonomic ranks and genome statistics, sorted by total abundance; 
               (2) a MultiQC-compatible report file for heatmap visualization 
               of relative bin abundances; (3) a filtered stats table for 
               MultiQC general table; (4) a JSON file for MultiQC scatter plot 
               of bin quality (completeness vs. contamination) with MIMAG-based 
               color coding. Unmapped reads are tracked as a separate 'unmapped_to_bin' 
               entry. Bins are extracted from individual FASTA files in a 
               dereplicated bins directory.
  Input files: 
    - List of samtools coverage TSV files (one per sample, with .coverage.tsv suffix)
    - List of samtools flagstat files (one per sample, with .flagstat suffix)
    - Directory of bin FASTA files (one .fasta per bin, named as bin_id.fasta)
    - DRep genomeInformation.csv file with completeness, contamination, length, N50
    - GTDB-Tk batch output file (user_genome, classification columns)
  Output files:
    - Main abundance table (TSV): per-bin metrics + taxonomy + per-sample abundances
    - MultiQC report file (TSV): relative abundances (%) for heatmap
    - MultiQC table file (TSV): filtered high-quality bins (completeness>50%, contamination<10%)
    - MultiQC JSON file: CheckM2 quality scatter plot data with MIMAG color categories
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

import pandas as pd
import numpy as np
import argparse
import json
import re
import os

################################################
# Function
################################################


def bins_contigs_compositions(folder):
    '''
    - Function: Associate each contig with the bin were the contig is retrieved.
    - Input: Bins folder (e.g. dereplicated bins folder)
    - Output: Dictionnary {contig_name : bin_1 }
    '''
    contigs_to_bins = dict()
    list_bins = list()
    for fi in os.listdir(folder):
        bin_fasta_file = os.path.join(folder, fi)
        with open(bin_fasta_file) as bin_file:
            bin_name = fi.rsplit('.', maxsplit=1)[0]
            if bin_name not in list_bins:
                list_bins.append(bin_name)
            for line in bin_file:
                if line.startswith('>'):
                    contig_name = line.split()[0].strip().lstrip('>')
                    contigs_to_bins[contig_name] = bin_name

    return contigs_to_bins, list_bins


def calculate_sample_bins_abundances(list_of_coverage_files, list_of_flagstat_files, contigs_to_bins, list_bins):
    '''
    - Function: Create a panda datadrame of bins abundances with sample as columns, and bins as rows.
                Adds a mean depth column that takes into account contigs lengths of bins.
                Flagstats files allows to add a row of proportion of unassigned reads.
                Finally, Adds a bin abundances sum column.
    - Input: -  list_of_coverage_files : samtools coverage files for each sample.
             -  list_of_flagstat_files : samtools flagstats files for each sample.
             -  contigs_to_bins : Dictionnary that link each contig with the bin were the contig is retrieve. \
                Create with bins_contigs_compositions() function previously.
             -  list_bins : Path of bins folder with fasta files.
    - Output: sample_to_bins_abundances :  Pandas dataframe with bins abundances per sample.
    '''
    list_samples = list()
    for coverage_file in list_of_coverage_files:
        sample_name = coverage_file.replace('.coverage.tsv', '')
        list_samples.append(sample_name)

    list_of_df_by_bin = list()
    for coverage_file in list_of_coverage_files:
        sample_name = coverage_file.replace('.coverage.tsv', '')
        df = pd.read_csv(coverage_file, sep='\t')

        df['bin'] = df['#rname'].apply(lambda c: contigs_to_bins[c])
        df['total_depth'] = df['meandepth'] * df['endpos']

        df_by_bin = df.groupby(['bin']).agg({"#rname": 'count',  #"#rname":';'.join,
                               'numreads': sum,
                               'endpos': sum,
                                "total_depth": sum,
                                      }).reset_index()

        df_by_bin = df_by_bin.rename(columns={"#rname": "contig_count", "endpos": 'bin_size'})
        df_by_bin = df_by_bin.set_index('bin')

        df_by_bin['meandepth'] = df_by_bin["total_depth"]/df_by_bin["bin_size"]
        df_by_bin = df_by_bin.drop(["total_depth"], axis=1)
        df_by_bin = df_by_bin.rename(columns={"numreads": f"numreads_{sample_name}", "meandepth": f"meandepth_{sample_name}"})
        list_of_df_by_bin.append(df_by_bin)
    sample_to_bins_abundances = pd.concat(list_of_df_by_bin, axis=1)
    # Add number of unassigned reads in the table
    flagstat_regexes = {
        "primary": r"(\d+) \+ (\d+) primary",
        "primary_mapped": r"(\d+) \+ (\d+) primary mapped \((.+):(.+)\)"
        }
    for flagstat_file in list_of_flagstat_files:
        sample_name = flagstat_file.replace('.flagstat', '')
        total = 0
        mapped = 0
        with open(flagstat_file, 'r') as fi:
            fi = "".join(fi.readlines())
            for flag, reg in flagstat_regexes.items():
                r_search = re.search(reg, fi)
                if r_search and flag == "primary":
                    total = int(r_search.group(0).split()[0])
                elif r_search and flag == "primary_mapped":
                    mapped = int(r_search.group(0).split()[0])

            sample_to_bins_abundances.loc['unmapped_to_bin', f"numreads_{sample_name}"] = \
            total - mapped

    columns_numreads = [c for c in sample_to_bins_abundances.columns if c.startswith('numreads')]
    columns_meandepth = [c for c in sample_to_bins_abundances.columns if c.startswith('meandepth')]
    sample_to_bins_abundances['sum_numreads'] = sample_to_bins_abundances[columns_numreads].sum(axis=1)
    sample_to_bins_abundances['sum_meandepth'] = sample_to_bins_abundances[columns_meandepth].sum(axis=1)
    sample_to_bins_abundances = sample_to_bins_abundances.drop(['bin_size'], axis=1)
    sample_to_bins_abundances = sample_to_bins_abundances.loc[:, ~sample_to_bins_abundances.columns.duplicated()]

    return sample_to_bins_abundances


def add_genomes_informations(genome_info, final_bins):
    '''
    - Function: Adds bins informations metrics (completeness,contamination,length,N50)
    - Input: drep/data_tables/genomeInformation.csv file generated previously.
    '''
    df_drep = pd.read_csv(genome_info)
    df_drep = df_drep.rename(columns = {'N50': 'genome_N50', 'length': 'genome_length'})
    df_drep['genome'] = df_drep['genome'].str.split('.').str[:-1].str.join('.')
    df_drep = df_drep.loc[df_drep['genome'].isin(final_bins)]
    df_drep = df_drep.drop(['centrality'], axis=1)
    df_drep = df_drep.set_index('genome')

    return df_drep


def return_lowest_taxo_rank(taxo):
    '''
    - Function: Returns the lowest non-null taxonomic affiliation from the entire taxonomic \
                obtained with gtdb-tk.
    - Input: d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium beijerinckii
    - Output: s__Clostridium beijerinckii
    '''
    ranks = ["s__", "g__", "f__", "o__", "c__", "p__", "d__"]
    if not type(taxo) == float:
            for i in range(-1, -len(ranks)-1, -1):
                if taxo.split(';')[i] != ranks[-i-1]:
                    return taxo.split(';')[i]
    return "unknown"


def add_bins_affiliations(affiliations_predictions_file):
    '''
    - Function: Reads gtdb-tk taxonomic affiliations file and add sample_name associated.
                Splits the classificiations column into differents rank taxonomic columns.
    '''
    affiliations = pd.read_csv(affiliations_predictions_file, sep='\t', usecols = ['user_genome', 'classification']).set_index('user_genome')
    affiliations['genome_name'] = affiliations['classification'].apply(lambda taxo: return_lowest_taxo_rank(taxo))
    ranks = {0: "Domain", 1: "Phylum", 2: "Class", 3: "Order", 4: "Family", 5: "Genus", 6: "Species"}
    for i, cur_rank in ranks.items():
        affiliations[cur_rank] = affiliations['classification'].str.split(';').str[i].str.split('__').str[-1]
    affiliations = affiliations[['genome_name'] + ['Domain'] + ["Phylum"] + ["Class"] + ["Order"] + ["Family"] + ["Genus"] + ["Species"]]

    return affiliations


def write_general_output_file(affiliations, informations, abundances, output_file):
    '''
    - Function: Concatanate abundances, genomes informations and affiliations DataFrames \
      made previously into one global table output.
      Sorts the bins by abundances.
    '''
    bins_general = pd.concat([affiliations, informations, abundances], axis=1)
    bins_general = \
    bins_general.sort_values('sum_numreads', ascending=False)
    bins_general['genome_id'] = bins_general.index
    bins_general = bins_general[['genome_id'] + [c for c in bins_general.columns if c != "genome_id"]]
    bins_general.loc[bins_general['genome_id'] == "unmapped_to_bin", 'genome_name'] = 'unmapped_to_bin'
    bins_general.to_csv(output_file, sep='\t', index=False)

    return bins_general


def write_report_file(general_table, report_file, checkm2_file, table_file):
    '''
    - Function: Write mqc files in order to make MultiQC output figures.
    '''
    report_df = general_table.copy()
    table_df = general_table.copy()
    idx = report_df.index.tolist()
    idx.remove('unmapped_to_bin')
    report_df = report_df.reindex(idx + ['unmapped_to_bin'])
    report_df = report_df.set_index('genome_name')
    report_cols = [col for col in report_df.columns if col.startswith('numreads') or col == "sum_numreads"]
    table_cols = [col for col in table_df.columns if not col.startswith('meandepth_') and not col.startswith('numreads_') \
    and not col == "Domain" and not col == "Phylum" and not col == "Class" and not col == "Order" \
    and not col == "Family" and not col == "Genus" and not col == "Species"]
    table_df = table_df[table_cols]
    table_df = table_df[(table_df['completeness'] > 50) & (table_df['contamination'] < 10)]
    report_df = report_df[report_cols]
    # Normalize library size by transforming values as samples percentages abundances
    for column in report_df.columns:
        report_df[column] = \
        report_df[column] / report_df[column].sum(axis=0) * 100
    report_df = report_df.transpose()
    ##
    # Filter heatmap to keeponly abundances genomes. We keep a genome if it's proportion is more than 5% in at least one sample
    # filter = ((report_df>=5).any()) | (report_df.columns == "unmapped_to_bin")
    unmapped_col = report_df["unmapped_to_bin"]
    report_df = report_df.iloc[:, :30]
    report_df.loc[:, "unmapped_to_bin"] = unmapped_col
    ##
    report_df.index.name = "sample"
    report_df = report_df.reset_index(level='sample')
    report_df.to_csv(report_file, sep='\t', index=False)
    table_df.to_csv(table_file, sep="\t", index=False)
    ### generate .json file for checkm2 quality bins scatterplot
    checkm_to_json = dict(id = 'bins_quality',
                    section_name = 'Bins Quality overview',
                    description = "Quality of bins in terms of completeness and contamination calculated by Checkm2. The points are colored according to their quality, according to the MIMAG standards defined previously (see Bins Counts quality section). Genomes with the best quality (100\% completeness and 0\% contamination) are located in the lower right corner of the graph. ",
                    plot_type = 'scatter',
                    anchor = 'bins_quality',
                    pconfig = dict(
                        title = 'Bins quality overview',
                        ylab = 'Contamination',
                        xlab = 'Completeness'))
    conditions = [
        ((general_table['completeness'] > 90) & (general_table['contamination'] < 5)),
        ((general_table['completeness'] > 50) & (general_table['contamination'] < 10)),
        ((general_table['completeness'] < 50) | (general_table['contamination'] > 10))
        ]
    # create a list of color according to the quality of the bins
    values = ['#D5ECC2', '#FFD3B4', '#ffd92f']
    general_table['color'] = np.select(conditions, values)

    checkm_to_json['data'] = dict()
    for index, row in general_table.iterrows():
        checkm_to_json['data'][row['genome_name']] = dict(x = row["completeness"], y = row["contamination"], color = row["color"])

    json_data = json.dumps(checkm_to_json)
    f = open(checkm2_file, 'w')
    f.write(json_data)
    ###


def parse_arguments():
    # Manage parameters.
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--list_of_coverage_files", nargs='+', action='append',\
    required = True, help = "List of samtools coverages file")

    parser.add_argument("-f", "--list_of_flagstats_files", nargs='+', action='append',\
    required = True, help = "List of samtools flagstats file")

    parser.add_argument('-b', '--bins_folder', required = True, help = \
    'fasta file containing sequences of contigs.')

    parser.add_argument('-g', '--genomes_informations', required = True, help = \
    'drep output file containing genomes informations.')

    parser.add_argument('-a', '--affiliations_predictions', required = True, help = \
    'gtdbtk taxonomic affiliations predictions output file.')

    parser.add_argument('-o', '--output_file', required = True, \
    help = 'Name of output file containing bins abundances per sample.')

    parser.add_argument('-r', '--report_file', required = True, \
    help = 'Name of report file as input for generate multiqc heatmap.')

    parser.add_argument('-c', '--checkm_file', required = True, \
    help = 'Name of checkm2 stats file as input for generate multiqc scatterplot.')

    parser.add_argument('-t', '--table_file', required = True, \
    help = 'Name of binning stats table file as input for generate multiqc table.')

    args = parser.parse_args()
    return args

###################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################


def main():

    args = parse_arguments()

    coverages_files = args.list_of_coverage_files[0]
    flagstats_files = args.list_of_flagstats_files[0]

    bins_compositions, list_bins = bins_contigs_compositions(args.bins_folder)

    bins_abundances = calculate_sample_bins_abundances(coverages_files, flagstats_files, bins_compositions, list_bins)

    bins_taxo_affiliations = add_bins_affiliations(args.affiliations_predictions)

    bins_informations = add_genomes_informations(args.genomes_informations, list_bins)

    general_table = write_general_output_file(bins_taxo_affiliations, bins_informations, bins_abundances, args.output_file)

    write_report_file(general_table, args.report_file, args.checkm_file, args.table_file)

if __name__ == '__main__':
    main()