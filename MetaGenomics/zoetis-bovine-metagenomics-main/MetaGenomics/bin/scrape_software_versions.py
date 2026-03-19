#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re
import os 

regexes = {
    'metagWGS': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'BWA-MEM2': ['v_bwa.txt', r"(\S+)"],
    'Minimap2': ['v_minimap.txt', r"(\S+)"],
    'Cutadapt': ['v_cutadapt.txt', r"(\S+)"],
    'Sickle': ['v_sickle.txt', r"sickle version (\S+)"],
    'KronaTools': ['v_kronatools.txt', r"KronaTools (\S+)"],
    'Python': ['v_python.txt', r"Python (\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'SPAdes': ['v_spades.txt', r"SPAdes genome assembler v(\S+)"],
    'Kaiju': ['v_kaiju.txt', r"Kaiju (\S+)"],
    'Kraken': ['v_kraken.txt', r"Kraken (\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],

}
results = OrderedDict()
results['metagWGS'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['Python'] = '<span style="color:#999999;\">N/A</span>'
results['Cutadapt'] = '<span style="color:#999999;\">N/A</span>'
results['Sickle'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['Kaiju'] = '<span style="color:#999999;\">N/A</span>'
results['Kraken'] = '<span style="color:#999999;\">N/A</span>'
results['KronaTools'] = '<span style="color:#999999;\">N/A</span>'
results['SPAdes'] = '<span style="color:#999999;\">N/A</span>'
results['BWA-MEM2'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'

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
section_href: 'https://lizard.bio'
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
