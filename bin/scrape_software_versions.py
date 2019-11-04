#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'nf-core/bcellmagic': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"(\S+)"],
    'vsearch': ['v_vsearch.txt', r"(\S+)"],
    'cd-hit': ['v_cdhit.txt', r"cd-hit (\S+)"],
    'muscle': ['v_muscle.txt', r"MUSCLE (\S+)"],
    'igblast': ['v_igblast.txt', r"igblast  (\S+)"],
    'airr': ['v_airr.txt', r"airr (\S+)"],
    'presto': ['v_presto.txt', r"(\S+)"],
    'changeo': ['v_changeo.txt', r"(\S+)"],
    'R': ['v_R.txt', r"(\S+)"],
    'r-shazam': ['v_shazam.txt', r"(\S+)"],
    'r-alakazam': ['v_alakazam.txt', r"(\S+)"],
    'r-tigger': ['v_tigger.txt', r"(\S+)"]
}
results = OrderedDict()
results['nf-core/bcellmagic'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Remove software set to false in results
for k in results:
    if not results[k]:
        del(results[k])

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'nf-core/bcellmagic Software Versions'
section_href: 'https://github.com/nf-core/bcellmagic'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
    for k,v in results.items():
        f.write("{}\t{}\n".format(k,v))
