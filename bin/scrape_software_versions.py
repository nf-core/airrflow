#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'nf-core/bcellmagic': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC (\S+)"],
    'MultiQC': ['v_multiqc.txt', r"version (\S+)"],
    'vsearch': ['v_vsearch.txt', r"vsearch v(\S+)"],
    'muscle': ['v_muscle.txt', r"MUSCLE v(\S+)"],
    'presto': ['v_presto.txt', r"(\S+)"],
    'changeo': ['v_changeo.txt', r"(\S+)"],
    'R': ['v_R.txt', r"version (\S+)"],
    'r-shazam': ['v_shazam.txt', r"(\S+)"],
    'r-alakazam': ['v_alakazam.txt', r"(\S+)"],
    'r-tigger': ['v_tigger.txt', r"(\S+)"]
}
results = OrderedDict()
results['nf-core/bcellmagic'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results["FastQC"] = '<span style="color:#999999;">N/A</span>'
results["MultiQC"] = '<span style="color:#999999;">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del results[k]

# Dump to YAML
print(
    """
id: 'software_versions'
section_name: 'nf-core/bcellmagic Software Versions'
section_href: 'https://github.com/nf-core/bcellmagic'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
"""
)
for k, v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out regexes as csv file:
with open("software_versions.csv", "w") as f:
    for k, v in results.items():
        f.write("{}\t{}\n".format(k, v))
