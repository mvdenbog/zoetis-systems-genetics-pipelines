#!/usr/bin/env python

"""----------------------------------------------------------------------------
  Script Name: collect_software_versions.py
  Description: Generates MultiQC-compatible software version reports by parsing 
               version output files from a metagenomic analysis pipeline. For 
               each supported tool (e.g., BWA-MEM2, Minimap2, SPAdes, CheckM2, 
               GTDB-Tk), the script: (1) checks for the existence of a version 
               file (e.g., v_bwa.txt); (2) applies a tool-specific regex pattern 
               to extract the version string; (3) formats the result as HTML 
               with version prefix. Tools without detected version files remain 
               marked as 'N/A' and are excluded from the final output. Outputs: 
               (1) YAML-formatted MultiQC custom content module (printed to 
               stdout) with an HTML definition list of tool versions; (2) a 
               tab-separated CSV file (software_versions.csv) listing detected 
               tools and their versions. Designed for integration into MultiQC 
               reports to document pipeline reproducibility.
  Input files: 
    - Version text files (v_*.txt) containing tool version output, one per 
      supported software (e.g., v_bwa.txt, v_spades.txt, v_checkm2.txt)
  Output files:
    - Stdout: YAML/HTML MultiQC custom content module (id: software_versions)
    - software_versions.csv: TSV file with columns [tool_name, version]
  Created By:  Mathias Vandenbogaert
  Date:        2026-05-06
-------------------------------------------------------------------------------
"""

# Metadata
__author__ = 'Mathias Vandenbogaert'
__version__ = '0.1'
__email__ = 'mathias.vandenbogaert@lizard.bio'
__status__ = 'dev'

from __future__ import print_function
from collections import OrderedDict
import re
import os 

regexes = {
    'metagWGS': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'BWA-MEM2': ['v_bwa.txt', r"(\S+)"],
    'Minimap2': ['v_minimap2.txt', r"(\S+)"],
    'Cutadapt': ['v_cutadapt.txt', r"(\S+)"],
    'Sickle': ['v_sickle.txt', r"sickle version (\S+)"],
    'KronaTools': ['v_kronatools.txt', r"KronaTools (\S+)"],
    'Python': ['v_python.txt', r"Python (\S+)"],
    'CD-HIT': ['v_cdhit.txt', r"CD-HIT version (\S+)"],
    'FeatureCounts': ['v_featurecounts.txt', r"featureCounts v(\S+)"],
    'Diamond': ['v_diamond.txt', r"diamond v(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'Megahit': ['v_megahit.txt', r"MEGAHIT v(\S+)"],
    'SPAdes': ['v_spades.txt', r"SPAdes genome assembler v(\S+)"],
    'Hifiasm': ['v_hifiasm_meta.txt', r"ha base version: (\S+)"],
    'MetaFlye': ['v_metaflye.txt', r"(\S+)"],
    'metaMDBG': ['v_metamdbg.txt', r"Version: (\S+)"],
    'Quast': ['v_quast.txt', r"QUAST v(\S+)"],
    'Kaiju': ['v_kaiju.txt', r"Kaiju (\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'Eggnog-Mapper': ['v_eggnogmapper.txt', r"emapper-(\S+)"],
    'Concoct': ['v_concoct.txt', r"concoct (\S+)"],
    'Maxbin': ['v_maxbin.txt', r"MaxBin (\S+)"],
    'Metabat2': ['v_metabat2.txt', r"version (\S+)"],
    'CheckM2': ['v_checkm2.txt', r"(\S+)"],
    'Binette': ['v_binette.txt', r"(\S+)"],
    'GTDBTK': ['v_gtdbtk.txt', r"...::: GTDB-Tk v(\S+)"],
    'dRep': ['v_dRep.txt', r"v(\S+)"],
    'tRNAscan-SE': ['v_tRNAscan.txt', r"tRNAscan-SE (\S+)"],
    'Barrnap': ['v_barrnap.txt', r"barrnap (\S+)"],
    'Prodigal': ['v_prodigal.txt', r"Prodigal V(\S+):"],
}
results = OrderedDict()
results['metagWGS'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['Python'] = '<span style="color:#999999;\">N/A</span>'
results['Cutadapt'] = '<span style="color:#999999;\">N/A</span>'
results['Sickle'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['Kaiju'] = '<span style="color:#999999;\">N/A</span>'
results['KronaTools'] = '<span style="color:#999999;\">N/A</span>'
results['Megahit'] = '<span style="color:#999999;\">N/A</span>'
results['SPAdes'] = '<span style="color:#999999;\">N/A</span>'
results['MetaFlye'] = '<span style="color:#999999;\">N/A</span>'
results['Hifiasm'] = '<span style="color:#999999;\">N/A</span>'
results['metaMDBG'] = '<span style="color:#999999;\">N/A</span>'
results['Quast'] = '<span style="color:#999999;\">N/A</span>'
results['BWA-MEM2'] = '<span style="color:#999999;\">N/A</span>'
results['Minimap2'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
results['Prokka'] = '<span style="color:#999999;\">N/A</span>'
results['Diamond'] = '<span style="color:#999999;\">N/A</span>'
results['CD-HIT'] = '<span style="color:#999999;\">N/A</span>'
results['FeatureCounts'] = '<span style="color:#999999;\">N/A</span>'
results['Eggnog-Mapper'] = '<span style="color:#999999;\">N/A</span>'
results['Concoct'] = '<span style="color:#999999;\">N/A</span>'
results['Metabat2'] = '<span style="color:#999999;\">N/A</span>'
results['Maxbin'] = '<span style="color:#999999;\">N/A</span>'
results['CheckM2'] = '<span style="color:#999999;\">N/A</span>'
results['Binette'] = '<span style="color:#999999;\">N/A</span>'
results['dRep'] = '<span style="color:#999999;\">N/A</span>'
results['GTDBTK'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'
results['tRNAscan-SE'] = '<span style="color:#999999;\">N/A</span>'
results['Barrnap'] = '<span style="color:#999999;\">N/A</span>'
results['Prodigal'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    if os.path.exists(v[0]):
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))

# Remove unused softwares
results = { k:v for k,v in results.items() if results[k]!='<span style="color:#999999;\">N/A</span>'}

# Dump to YAML
print('''
id: 'software_versions'
section_name: 'metagWGS Software Versions'
section_href: 'https://forge.inrae.fr/genotoul-bioinfo/metagwgs  '
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k, v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
    for k, v in results.items():
        f.write("{}\t{}\n".format(k, v))