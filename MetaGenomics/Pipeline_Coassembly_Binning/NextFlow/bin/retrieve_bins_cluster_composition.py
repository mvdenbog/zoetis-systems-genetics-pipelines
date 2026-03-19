#!/usr/bin/env python

"""----------------------------------------------------------------------------
  Script Name: drep_cluster_to_bins.py
  Description: Processes dRep dereplication output to create a mapping table 
               linking each representative bin (genome) to all member bins in 
               its cluster. Reads two dRep data_tables CSV files: (1) Wdb.csv 
               containing cluster representatives (genome, cluster columns); 
               (2) Cdb.csv containing cluster compositions (genome, 
               secondary_cluster columns). For each secondary cluster, 
               aggregates all member genome IDs into a comma-separated list, 
               then joins with the representative mapping to produce a two-
               column TSV: representative bin ID and comma-delimited list of 
               all bins (including representative) belonging to that cluster. 
               Useful for tracking dereplication relationships, aggregating 
               metrics at the representative level, or filtering analyses to 
               primary clusters only.
  Input files: 
    - Wdb.csv: dRep representatives file (columns: genome, cluster)
    - Cdb.csv: dRep cluster composition file (columns: genome, secondary_cluster)
  Output files:
    - TSV table with two columns: 'representative' (cluster representative 
      genome ID) and 'composition' (comma-separated list of all member genome 
      IDs in that cluster)
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
import sys, re
import argparse
import pandas as pd

def read_representatives_file(fi):
    df = pd.read_csv(fi)
    df = df[["genome","cluster"]]
    df = df.set_index('cluster')
    return df

def read_clusters_file(fi, df_repr, output_fi):
    df_clust = pd.read_csv(fi)
    df_clust = df_clust[["genome","secondary_cluster"]]
    composition = df_clust.groupby(['secondary_cluster'])['genome'].apply(list)
    df_clust = composition.to_frame()['genome'].apply(lambda x: ",".join(x) ).to_frame()
    df = pd.merge(df_clust, df_repr,  left_index=True, right_index=True)
    df = df.rename(columns={"genome_y": "representative", "genome_x": "composition"})
    df = df[["representative","composition"]]

    df.to_csv(output_fi ,sep='\t', index=False)

def parse_arguments():
    # Manage parameters.
    parser = argparse.ArgumentParser(description='Process dRep output to map representative bins to cluster members.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-c', '--clusters', required = True, help = \
    'data_tables/Cdb.csv dRep output file of bins cluster compositions.')

    parser.add_argument('-r', '--representatives', required = True, help = \
    'data_tables/Wdb.csv dRep output file of representative bins clusters.')

    parser.add_argument('-o', '--output_file', required = True, \
    help = 'representative-cluster_to_bins output file.')

    args = parser.parse_args()
    return args

###################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################


def main():

    args = parse_arguments()

    df_repr = read_representatives_file(args.representatives)

    read_clusters_file(args.clusters, df_repr, args.output_file)

if __name__ == '__main__':
    main()