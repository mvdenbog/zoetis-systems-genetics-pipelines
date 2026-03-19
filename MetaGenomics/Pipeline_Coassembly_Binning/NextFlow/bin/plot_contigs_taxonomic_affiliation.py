#!/usr/bin/env python3

"""----------------------------------------------------------------------------
  Script Name: plot_contigs_taxonomic_affiliation.py
  Description: Generates interactive Plotly bar plots visualizing taxonomic 
               composition of metagenomic samples based on contig-level 
               taxonomic assignments. Input files are per-sample quantification 
               tables (from merge_contig_quantif_perlineage.py) containing 
               contig abundances (read counts, mean depth) and taxonomic 
               lineages. Produces: (1) stacked bar plots showing relative 
               abundance (% reads, % depth) of taxonomic ranks (Kingdom to 
               Species) per sample; (2) interactive tabbed HTML report with 
               per-rank bar plots of the top-N most abundant taxa (default: 10), 
               with "Unknown" and "Other <rank>" categories for unassigned or 
               low-abundance taxa. Taxa are colored using a fixed palette; 
               "Unknown" and "Other" categories use distinct colors. Lineage 
               strings are parsed to determine the lowest resolved taxonomic 
               rank for grouping. Supports multiple samples with consistent 
               color mapping and sample ordering. Output is a directory of 
               HTML files viewable in any modern browser.
  Input files: 
    - One or more <sample>_quantif_percontig.tsv files (from 
      merge_contig_quantif_perlineage.py) with columns: lineage_by_level, 
      tax_id_by_level, nb_reads, depth, nb_contigs, etc.
  Output files (in --output_dir, default: 'plots'):
    - abundance_per_rank.html: stacked bar plots of rank-level abundances 
      (% reads and % depth) across samples
    - most_abundant_taxa.html: interactive tabbed HTML with per-rank bar plots 
      of top-N taxa per sample (tabs: Kingdom through Species)
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
import sys
import plotly.express as px
import pandas as pd
import os
from collections import defaultdict

def get_rank(lineage_str, ranks):
    
    unknown_lineage =['None; None; None; None; None; None; uncultured prokaryote',
                     ' Unable to find taxonomy consensus',
                     'No affiliation',
                     'None; None; None; None; None; None; uncultured organism']

    if lineage_str in unknown_lineage:
        return "Unknown"
    
    lineage = lineage_str.split(";")
    lineage_len = len(lineage)
    for t in lineage[::-1]:
        if t.strip() == "None":
            lineage_len -= 1 
        else:
            break
    return ranks[lineage_len - 1 ]


def get_top_taxa(df, taxon_col, abd, n, ignored_taxon):
    df_rank_grouped = df.groupby([taxon_col]).agg({
                                        "nb_reads":"sum", 
                                        "depth":"mean", 
                                        abd:"sum", }).reset_index()
    df_rank_grouped =  df_rank_grouped.sort_values(by=abd, ascending=False)
    return list(df_rank_grouped.loc[df_rank_grouped[taxon_col] != ignored_taxon][taxon_col].head(n))


def write_tab_html(html_figs, outdir):
    html_head = """<!DOCTYPE html>
<html>
<head>
<meta name="viewport" content="width=device-width, initial-scale=1">

<style>
body {font-family: Arial;}

/* Style the tab */
.tab {
  overflow: hidden;
  border: 1px solid #ccc;
  background-color: #f1f1f1;
}

/* Style the buttons inside the tab */
.tab button {
  background-color: inherit;
  float: left;
  border: none;
  outline: none;
  cursor: pointer;
  padding: 14px 16px;
  transition: 0.3s;
  font-size: 17px;
}

/* Change background color of buttons on hover */
.tab button:hover {
  background-color: #ddd;
}

/* Create an active/current tablink class */
.tab button.active {
  background-color: #ccc;
}

/* Style the tab content */
.tabcontent {
  display: none;
  padding: 6px 12px;
  border: 1px solid #ccc;
  border-top: none;
}

</style>
</head>
<body>
"""


    html_tail = """
<script>
function openTab(evt, taxonName) {
  var i, tabcontent, tablinks;
  tabcontent = document.getElementsByClassName("tabcontent");
  for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
  }
  tablinks = document.getElementsByClassName("tablinks");
  for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }
  document.getElementById(taxonName).style.display = "block";
  evt.currentTarget.className += " active";
}

// Get the element with id="defaultOpen" and click on it
document.getElementById("defaultOpen").click();
</script>
   
</body>
</html> 
"""

    tab_div = """
<div class="tab">
    <button class="tablinks" onclick="openTab(event, 'Kingdom')" id="defaultOpen">Kingdom</button>
    <button class="tablinks" onclick="openTab(event, 'Phylum')">Phylum</button>
    <button class="tablinks" onclick="openTab(event, 'Class')">Class</button>
    <button class="tablinks" onclick="openTab(event, 'Order')">Order</button>
    <button class="tablinks" onclick="openTab(event, 'Family')">Family</button>
    <button class="tablinks" onclick="openTab(event, 'Genus')">Genus</button>
    <button class="tablinks" onclick="openTab(event, 'Species')">Species</button>
  </div>
"""
    # for abd, html_figs_list in html_figs.items():
    #     abd_cleaned = abd.replace(' ', '_').replace('%', '').replace('_(', '').replace(')', '')
    #     outfile_name = f"plot_most_{abd_cleaned}_taxa.html"
    #     outfile = os.path.join(outdir, outfile_name)
    #     with open(outfile, "w") as fl:
    #         fl.write(html_head)
    #         fl.write(tab_div)
    #         for rank, html_fig in html_figs_list:
    #             fl.write(f'<div id="{rank}" class="tabcontent">\n')
    #             fl.write(html_fig)
    #             fl.write(f'</div>\n')
                
    #         fl.write(html_tail)

    rank_to_html_figs = defaultdict(list)
    for abd, html_figs_list in html_figs.items():
        for rank, html_fig in html_figs_list:
            rank_to_html_figs[rank].append(html_fig)

    outfile_name = f"most_abundant_taxa.html"
    outfile = os.path.join(outdir, outfile_name)

    with open(outfile, "w") as fl:
        fl.write(html_head)
        fl.write(tab_div)
        for rank, html_figs in rank_to_html_figs.items():
            fl.write(f'<div id="{rank}" class="tabcontent">\n')
            for html_fig in html_figs:
                fl.write(html_fig)
            fl.write(f'</div>\n')
            
        fl.write(html_tail)

def parse_affi_files(contig_affi_files):
    df_list = []
    for affi_file in contig_affi_files:
        df = pd.read_csv(affi_file, sep='\t')
        df['Abundance (% reads)'] = 100 * df['nb_reads'] / df['nb_reads'].sum()
        
        df['Abundance (% depth)'] = 100 * df['depth'] / df['depth'].sum()

        sample = os.path.basename(affi_file.replace('_quantif_percontig.tsv', ""))
        df['sample'] = sample
        df_list.append(df)
        
    df = pd.concat(df_list)
    return df


def make_rank_plot(df, abundance_types, rank2color, ranks, samples_sorted, outdir):

    df_rank = df.groupby(['rank', "sample"]).agg({"nb_reads":sum, 
                                                    "depth":"mean",
                                                    'Abundance (% depth)':"sum", 
                                             "Abundance (% reads)":"sum"}).reset_index()



    df_rank['abd'] = df_rank["Abundance (% reads)"].round(3).astype(str) + '%'

    fig = px.bar(df_rank, x="sample", y="Abundance (% reads)", color="rank",
                color_discrete_map=rank2color, template="plotly_white",
                category_orders={"rank":ranks, "sample":samples_sorted}, text="abd", 
                title='Abundance in read percentage of taxonomic ranks computed from contigs affiliation')
    html_fig_reads = fig.to_html(full_html=False, include_plotlyjs=True)
    


    df_rank['abd'] = df_rank['Abundance (% depth)'].round(3).astype(str) + '%'

    fig = px.bar(df_rank, x="sample", y='Abundance (% depth)', color="rank",
                color_discrete_map=rank2color, template="plotly_white",
                category_orders={"rank":ranks, "sample":samples_sorted}, text="abd", 
                title='Abundance in depth percentage of taxonomic ranks computed from contigs affiliation')

    html_fig_depth = fig.to_html(full_html=False, include_plotlyjs=False)

    outfile = os.path.join(outdir, f"abundance_per_rank.html")

    with open(outfile, "w") as fl:
        fl.write(html_fig_reads)      
        fl.write(html_fig_depth)
        
    # fig.write_html(f"Abundance_per_rank_{abd_cleaned}_cdn.html", include_plotlyjs='cdn')

def makes_abundant_taxa_plots(df, abd, ranks, samples_sorted, unknown_color, other_color, color_pal, top_n_taxon):
    
    rank2fig = {}
    for rank in ranks:

        rank_name = f'{rank}_name'
        rank_taxid = f'{rank}_taxid'
        df_to_grp = df.sort_values(by=abd, ascending=False)
  
        unkown_taxa = [" None", "None", ' Unable to find taxonomy consensus', None]
        df_to_grp.loc[df_to_grp[rank_name].isin(unkown_taxa), rank_name] = f"Unknown"

        df_to_grp.loc[df_to_grp[rank_name] == "Unknown", rank_taxid] = f"Unknown"

        
        top_taxa = get_top_taxa(df_to_grp, rank_name, abd, top_n_taxon, f"Unknown") 

        top_taxa.append(f"Unknown")

        df_to_grp = df_to_grp.sort_values(by=abd, ascending=False)


        df_to_grp.loc[~df_to_grp[rank_name].isin(top_taxa), rank_taxid] = f"Other {rank}"
        df_to_grp.loc[~df_to_grp[rank_name].isin(top_taxa), rank_name] = f"Other {rank}"

        

        df_to_grp = df_to_grp.sort_values(by=abd, ascending=False)
        

        df_rank_grouped = df_to_grp.groupby([rank_name, rank_taxid, 'sample',
                                        ]).agg({
                                                "nb_reads":"sum", 
                                                "depth":"mean", 
                                                abd:"sum", }).reset_index()

        if abd == 'Abundance (% reads)':
            title = f"Abundance in reads percent of the {top_n_taxon} most abundant {rank} computed from contigs affiliation"
        elif abd == 'Abundance (% depth)':
            title = f"Abundance in depth percent of the {top_n_taxon} most abundant {rank} computed from contigs affiliation"
        else:
            title = f"Abundance of the {top_n_taxon} most abundant {rank} computed from contigs affiliation"

        df_rank_grouped['abd'] = df_rank_grouped[abd].round(2).astype(str) + '%'

        df_rank_grouped['abundance_to_sort'] = df_rank_grouped[abd]
        df_rank_grouped.loc[df_rank_grouped[rank_name] == f"Unknown", 'abundance_to_sort'] = -2
        df_rank_grouped.loc[df_rank_grouped[rank_name] == f"Other {rank}", 'abundance_to_sort'] = -1

        df_rank_grouped = df_rank_grouped.sort_values(by='abundance_to_sort', ascending=False)
        
        df_rank_grouped = df_rank_grouped.rename(columns={rank_name: rank})

        fig = px.bar(df_rank_grouped, x="sample", y=abd,
                    color=rank, category_orders={"sample":samples_sorted},template="plotly_white",
                    color_discrete_map={f"Unknown":unknown_color, 
                                f"Other {rank}":other_color},
                    color_discrete_sequence=color_pal, text="abd", title=title)

        # fig.update_layout(
        #     width=210 * len(samples_sorted),
        #     height=600,)

        rank2fig[rank] = fig

    return rank2fig


def parse_arguments():
    """Parse script arguments."""
    parser = ArgumentParser(description="Generate interactive taxonomic composition plots from contig affiliation tables.",
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('affi_taxo_quantif', nargs='+', help='Taxonomic affiliation and quantitification file (<sample>_quantif_percontig.tsv) generated by merge_contig_quantif_perlineage.py script')

    parser.add_argument('--output_dir', default='plots', help="Name of the output directory")

    parser.add_argument('--nb_top_taxon', default=10, type=int, help="Plot only the top n most abundant taxa.")

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



    contig_affi_files = args.affi_taxo_quantif
    top_n_taxon = args.nb_top_taxon
    output_dir = args.output_dir

    os.makedirs(output_dir, exist_ok=True)

    ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

    no_affi_categories = ['No affiliation', 'Unable to find taxonomy consensus', 'Unknown']

    # COLOR PALETTE

    unknown_color="rgb(129,126,104)"

    rank2color = {r:c for r, c in  zip(ranks, px.colors.qualitative.Prism[1:])}
    rank2color['Unknown'] = unknown_color

    # import seaborn as sns
    # tab20 = list(sns.color_palette("tab20"))
    # tab20 = [f"rgb({r*250:.0f}, {g*250:.0f}, {b*250:.0f})" for r,g,b in tab20]
    color_pal = ['rgb(30, 117, 176)', 'rgb(171, 195, 227)', 'rgb(250, 125, 14)', 'rgb(250, 183, 118)', 
            'rgb(43, 157, 43)', 'rgb(149, 219, 135)', 'rgb(210, 38, 39)', 'rgb(250, 149, 147)', 
            'rgb(145, 101, 185)', 'rgb(193, 173, 209)', 'rgb(137, 84, 74)', 'rgb(192, 153, 145)', 
            'rgb(223, 117, 190)', 'rgb(242, 178, 206)', 'rgb(125, 125, 125)', 'rgb(195, 195, 195)', 
            'rgb(184, 185, 33)', 'rgb(215, 215, 138)', 'rgb(23, 186, 203)', 'rgb(155, 214, 225)']
    other_color = "rgb(190,191,195)"


    abundance_types = ['Abundance (% reads)', 'Abundance (% depth)']

    # PARSE contig_affi_files
    df = parse_affi_files(contig_affi_files)
    df['rank'] = df['lineage_by_level'].apply(lambda x: get_rank(x, ranks))

    samples_sorted = sorted(list(df['sample'].unique()))

    df.loc[df['lineage_by_level'].isin(no_affi_categories) ,"tax_id_by_level"] = "1; None; None; None; None; None; None"
    df.loc[df['lineage_by_level'].isin(no_affi_categories) ,"lineage_by_level"] = "Unknown; None; None; None; None; None; None"

    # Make ranks plots 
    make_rank_plot(df, abundance_types, rank2color, ranks, samples_sorted, output_dir)



    # Make n most abundant taxo plot

    ranks_names = [f"{r}_name" for r in ranks]
    ranks_taxid = [f"{r}_taxid" for r in ranks]

    df[ranks_taxid] = df["tax_id_by_level"].str.split(pat=";",expand=True)
    df[ranks_names] = df["lineage_by_level"].str.split(pat=";",expand=True)
    
    abd2html_figs = {}

    for abd in abundance_types:
        rank2fig = makes_abundant_taxa_plots(df,abd, ranks, 
                            samples_sorted, unknown_color, other_color, color_pal, top_n_taxon)
        html_figs = []
        for i, (rank, fig) in enumerate(rank2fig.items()):
            
            include_plotly = i == 0

            html_fig = fig.to_html(full_html=False, include_plotlyjs=include_plotly)
            html_figs.append((rank, html_fig )) 

        abd2html_figs[abd] = html_figs
        
        write_tab_html(abd2html_figs, output_dir)

    
if __name__ == '__main__':
    main()