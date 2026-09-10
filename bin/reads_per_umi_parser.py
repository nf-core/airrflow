#!/usr/bin/env python3
# Written by Jana Hoffmann and released under the MIT license (2026).

# reads_per_umi_parser.py
# Parse RNAseq and AIRRseq reads per UMI info

import pandas as pd
import glob
import argparse

parser = argparse.ArgumentParser(
    description="Get reads per UMI information for all samples."
)
parser.add_argument(
    "-rs",
    "--rnaseq",
    help="RNAseq workflow",
    action="store_true",
)
parser.add_argument(
    "-as",
    "--airrseq",
    help="AIRRseq workflow.",
    action="store_true",
)
args = parser.parse_args()

if args.rnaseq:
    sample_dfs = []
    for i, file in enumerate(glob.glob("*.csv")):
        sample_df = pd.read_csv(file, header=None, index_col=[1]).rename(columns={0:f"{file.split('.csv')[0]}"})
        sample_dfs.append(sample_df)
    reads_per_umi_df = pd.concat(sample_dfs, axis=1)

if args.airrseq:
    sample_counts = []
    for i, file in enumerate(glob.glob('UMI_clusterinfo/*R1_table.tab')):
        df = pd.read_csv(file, sep='\t')
        sample_counts.append(df['SEQCOUNT'].value_counts().rename(file.split('/')[-1].split('_R1_table.tab')[0]))
    reads_per_umi_df = pd.concat(sample_counts, axis=1)


reads_per_umi_df = reads_per_umi_df.reindex([i for i in range(reads_per_umi_df.index.min(),reads_per_umi_df.index.max()+1)])
reads_per_umi_df = reads_per_umi_df.fillna(0)
reads_per_umi_df.to_csv("reads_per_UMI.csv", index=True, index_label="")
