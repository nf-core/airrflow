#!/usr/bin/env python3
# Released under the MIT license (2026).

# log.py
# Parsing log files after each filtering step.

import glob
import os
import re
import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    description="Parse logs to identify the number of sequences passing through every step."
)
args = parser.parse_args()

# mapping of process to filename
processes = {
    "AddMetadata" : "_add-meta_command_log",
    "AssignGenes-igblast" : "_changeo_assigngenes_command_log",
    "CollapseDuplicates" : "_collapse_command_log",
    "FilterQuality" : "_fq_command_log",
    "FilterJunctionMod3" : "_jmod3_command_log",
    "MakeDB-igblast" : "_makedb_command_log",
    "ParseDb-split" : "_split_command_log",
    "SingleCellQC" : "_scqc_command_log",
    "CreateGermlines" : "_create-germlines_command_log",
    "RenameFile" : "_rename_command_log",
    "ConvertDb-fasta" : "_convertdb_command_log",
    "RemoveChimeric" : "_chimeric_command_log",
    "BulkOverlap" : "_contamination_command_log"
    #"ClonePass" : ""
}

columns = ["sample_id", "input", "input_size", "task", "output", "output_size"]
rows = []

# iterate over all log files
for log_file in sorted(glob.glob("logs/*command_log*")):
    with open(log_file) as f:
        # parse log file
        fields = dict(re.findall(r"([A-Z][A-Z0-9_]*)> *(\S*)", f.read()))
        task = fields["START"] + ("-" + fields["COMMAND"] if "COMMAND" in fields else "")

        sample_id = os.path.basename(log_file).split(processes[task])[0]
        if "BulkOverlap" in task:
            sample_id = None

        # get number of output files from process
        n_max = max([int(num) for key in fields.keys() for num in re.findall(r"[1-9]\d*", key)] + [0])

        if n_max == 0 or "ParseDb" in task:
            if "MakeDB" in task:
                input_file = fields['ALIGNER_FILE']
            else:
                input_file = fields.get("FILE")
            if "ParseDb" in task:
                input_size = fields["RECORDS"]
                # iterate all outputs (PASS1, PASS2, ...)
                for i in range(n_max):
                    output_size = fields["SIZE" + str(i + 1)]
                    output_file = fields["OUTPUT" + str(i + 1)]
                    rows.append([sample_id, input_file, input_size, task, output_file, output_size])
            else:
                output_size = fields["PASS"]
                output_file = fields["OUTPUT"]
                fail_size = (fields["FAIL"] if "FAIL" in fields else None)
                input_size = int(output_size) + (int(fail_size) if fail_size else 0)
                rows.append([sample_id, input_file, input_size, task, output_file, output_size])
        else:
            # iterate all outputs (PASS1, PASS2, ...)
            for i in range(n_max):
                output_size = fields["PASS" + str(i + 1)]
                fail_size = (fields["FAIL" + str(i + 1)] if "FAIL" + str(i + 1) in fields else None)
                input_size = int(output_size) + (int(fail_size) if fail_size else 0)
                input_file = fields["FILE" + str(i + 1)]
                rows.append([sample_id, input_file, input_size, task, fields["OUTPUT" + str(i + 1)], output_size])


# convert to dataframe and save as tsv
df = pd.DataFrame(rows, columns=columns)

# get sample id for SingleCellQC
df_log = df.copy()
df_log = df_log[df_log['task']=='SingleCellQC']
df_log['sample_id'] = df_log['input'].apply(lambda in_file: in_file.split('_meta-pass')[0])
df_log['task'] = 'HeavyChainFilter'

# add step for removed contaminants
if os.path.isfile('p-sc-duplicates.csv'):

    # get number of removed contaminants
    df_contamination = pd.read_csv('p-sc-duplicates.csv')
    before_contamination = df_contamination['sample_id'].value_counts()
    before_contamination.name = 'input_size'
    after_contamination = df_contamination[~df_contamination['sc_duplicate_cell']]['sample_id'].value_counts()
    after_contamination.name = 'output_size'
    df_contamination = pd.concat([before_contamination, after_contamination], axis=1).reset_index().rename(columns={'index':'sample_id'})
    df_contamination['task'] = 'RemoveContaminants'
    n_before_contamination = df_contamination.set_index('sample_id').to_dict()['input_size']

    # number of sequences before removing contaminants
    df_log['output_size'] = df_log['sample_id'].apply(lambda sample: n_before_contamination[sample] if sample in n_before_contamination.keys() else 0)
    df_combined = pd.concat([df_contamination, df_log, df[df['task']!='SingleCellQC']])
else:

    # update sample id
    df_combined = pd.concat([df_log, df[df['task']!='SingleCellQC']])

df_combined.to_csv("Table_sequences_assembled.tsv", sep="\t", index=False)
