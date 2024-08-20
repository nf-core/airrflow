#!/usr/bin/env python3
# Written by David Ladd, Gisela Gabernet and released under the MIT license (2020).

# log_parsing_no-umi.py
# Parsing log files for each of the steps for QC analysis.

import pandas as pd
import subprocess
import re

# Processes
processes = [
    "assemble_pairs",
    "filter_by_sequence_quality",
    "mask_primers",
    "deduplicates",
    "igblast",
]

# Path of logs will be:
# process_name/sample_name_command_log.txt

df_process_list = []

for process in processes:
    find = subprocess.check_output(["find", process, "-name", "*command_log*"])
    log_files = find.decode().split("\n")
    log_files = list(filter(None, log_files))

    if process in ["assemble_pairs"]:
        s_code = []
        pairs = []
        pass_pairs = []
        fail_pairs = []
        process_name = []

        for logfile in log_files:
            with open(logfile, "r") as f:
                for line in f:
                    if " START>" in line:
                        s_code.append(logfile.split("/")[1].split("_command_log")[0])
                        process_name.append(process)
                    elif "PAIRS>" in line:
                        pairs.append(line.strip().removeprefix("PAIRS> "))
                    elif "PASS>" in line:
                        pass_pairs.append(line.strip().removeprefix("PASS> "))
                    elif "FAIL>" in line:
                        fail_pairs.append(line.strip().removeprefix("FAIL> "))

        df_process = pd.DataFrame.from_dict(
            {
                "Sample": s_code,
                "start_pairs": pairs,
                "pass_pairs": pass_pairs,
                "fail_pairs": fail_pairs,
                "process": process_name,
            }
        )

        df_process_list.append(df_process)

    if process in ["filter_by_sequence_quality"]:
        s_code = []
        pass_pairs = []
        fail_pairs = []
        process_name = []

        for logfile in log_files:
            with open(logfile, "r") as f:
                for line in f:
                    if " START>" in line:
                        s_code.append(logfile.split("/")[1].split("_command_log")[0])
                        process_name.append(process)
                    elif "PASS>" in line:
                        pass_pairs.append(line.strip().removeprefix("PASS> "))
                    elif "FAIL>" in line:
                        fail_pairs.append(line.strip().removeprefix("FAIL> "))

        df_process = pd.DataFrame.from_dict(
            {
                "Sample": s_code,
                "pass_pairs": pass_pairs,
                "fail_pairs": fail_pairs,
                "process": process_name,
            }
        )

        df_process_list.append(df_process)

    elif process in ["mask_primers"]:
        s_code = []
        s_readtype = []
        output_file = []
        n_seqs = []
        n_pass = []
        n_fail = []
        process_name = []

        for logfile in log_files:
            if "_R1" in logfile:
                s_readtype.append("R1")
            elif "_R2" in logfile:
                s_readtype.append("R2")
            with open(logfile, "r") as f:
                for line in f:
                    if " START>" in line:
                        s_code.append(logfile.split("/")[1].split("_command_log")[0])
                        process_name.append(process)
                    elif "SEQUENCES>" in line:
                        n_seqs.append(line.strip().removeprefix("SEQUENCES> "))
                    elif "PASS>" in line:
                        n_pass.append(line.strip().removeprefix("PASS> "))
                    elif "FAIL>" in line:
                        n_fail.append(line.strip().removeprefix("FAIL> "))

        df_process = pd.DataFrame.from_dict(
            {
                "Sample": s_code,
                "readtype": s_readtype,
                "start": n_seqs,
                "pass": n_pass,
                "fail": n_fail,
                "process": process_name,
            }
        )

        df_process_list.append(df_process)

    elif process in ["pair_sequences", "repair_mates"]:
        s_code = []
        output_file = []
        seqs1 = []
        seqs2 = []
        pass_pairs = []
        process_name = []

        for logfile in log_files:
            with open(logfile, "r") as f:
                # print(f.read())
                for line in f:
                    if " START>" in line:
                        s_code.append(logfile.split("/")[1].split("_command_log")[0])
                        process_name.append(process)
                    elif "SEQUENCES1>" in line:
                        seqs1.append(line.strip().removeprefix("SEQUENCES1").removeprefix("> "))
                    elif "SEQUENCES2>" in line:
                        seqs2.append(line.strip().removeprefix("SEQUENCES2").removeprefix("> "))
                    elif "PASS>" in line:
                        pass_pairs.append(line.strip().removeprefix("PASS> "))

        df_process = pd.DataFrame.from_dict(
            {
                "Sample": s_code,
                "seqs_1": seqs1,
                "seqs_2": seqs2,
                "pass_pairs": pass_pairs,
                "process": process_name,
            }
        )

        df_process_list.append(df_process)

    elif process in ["deduplicates"]:
        s_code = []
        seqs = []
        unique = []
        duplicate = []
        undetermined = []
        process_name = []

        for logfile in log_files:
            with open(logfile, "r") as f:
                # print(f.read())
                for line in f:
                    if " START>" in line:
                        s_code.append(logfile.split("/")[1].split("_command_log")[0])
                        process_name.append(process)
                    elif "SEQUENCES>" in line:
                        seqs.append(line.strip().removeprefix("SEQUENCES> "))
                    elif "UNIQUE>" in line:
                        unique.append(line.strip().removeprefix("UNIQUE> "))
                    elif "DUPLICATE>" in line:
                        duplicate.append(line.strip().removeprefix("DUPLICATE> "))
                    elif "UNDETERMINED>" in line:
                        undetermined.append(line.strip().removeprefix("UNDETERMINED> "))

        df_process = pd.DataFrame.from_dict(
            {
                "Sample": s_code,
                "start_seqs": seqs,
                "unique": unique,
                "duplicate": duplicate,
                "undetermined": undetermined,
                "process": process_name,
            }
        )

        df_process_list.append(df_process)

    elif process in ["igblast"]:
        s_code = []
        pass_blast = []
        fail_blast = []
        for logfile in log_files:
            with open(logfile, "r") as f:
                # print(f.read())
                for line in f:
                    if "PASS>" in line:
                        s_code.append(logfile.split("/")[1].split("_command_log")[0])
                        pass_blast.append(line.strip().removeprefix("PASS> "))
                    elif "FAIL>" in line:
                        fail_blast.append(line.strip().removeprefix("FAIL> "))

        pass_fail = [list(map(int, pass_blast)), list(map(int, fail_blast))]
        repres_2 = [sum(x) for x in zip(*pass_fail)]

        df_process = pd.DataFrame.from_dict(
            {
                "Sample": s_code,
                "repres_2": repres_2,
                "pass_igblast": pass_blast,
                "fail_igblast": fail_blast,
            }
        )

        df_process_list.append(df_process)

    elif process in ["define_clones"]:
        s_code = []
        seqs = []
        clones = []
        pass_clones = []
        fail_clones = []
        process_name = []

        for logfile in log_files:
            with open(logfile, "r") as f:
                # print(f.read())
                for line in f:
                    if " START>" in line:
                        s_code.append(logfile.split("/")[1].split("_command_log")[0])
                        process_name.append(process)
                    elif "RECORDS>" in line:
                        seqs.append(line.strip().removeprefix("RECORDS> "))
                    elif "CLONES>" in line:
                        clones.append(line.strip().removeprefix("CLONES> "))
                    elif "PASS>" in line:
                        pass_clones.append(line.strip().removeprefix("PASS> "))
                    elif "FAIL>" in line:
                        fail_clones.append(line.strip().removeprefix("FAIL> "))

        df_process = pd.DataFrame.from_dict(
            {
                "Sample": s_code,
                "start_seqs": seqs,
                "clones": clones,
                "pass_clones": pass_clones,
                "fail_clones": fail_clones,
                "process": process_name,
            }
        )

        df_process_list.append(df_process)

    elif process in ["create_germlines"]:
        s_code = []
        seqs = []
        pass_clones = []
        fail_clones = []
        process_name = []

        for logfile in log_files:
            with open(logfile, "r") as f:
                # print(f.read())
                for line in f:
                    if " START>" in line:
                        s_code.append(logfile.split("/")[1].split("_command_log")[0])
                        process_name.append(process)
                    elif "RECORDS>" in line:
                        seqs.append(line.strip().removeprefix("RECORDS> "))
                    elif "PASS>" in line:
                        pass_clones.append(line.strip().removeprefix("PASS> "))
                    elif "FAIL>" in line:
                        fail_clones.append(line.strip().removeprefix("FAIL> "))

        df_process = pd.DataFrame.from_dict(
            {
                "Sample": s_code,
                "start_seqs": seqs,
                "pass_clones": pass_clones,
                "fail_clones": fail_clones,
                "process": process_name,
            }
        )

        df_process_list.append(df_process)

colnames = [
    "Sample",
    "Sequences",
    "Assemble_pairs",
    "Filtered_quality",
    "Mask_primers_R1",
    "Mask_primers_R2",
    "Unique",
    "Representative_2",
    "Igblast",
]

# Tables provide extra info and help debugging
df_process_list[0].to_csv(
    path_or_buf="Table_all_details_assemble_mates.tsv",
    sep="\t",
    header=True,
    index=False,
)
df_process_list[1].to_csv(
    path_or_buf="Table_all_details_filter_quality.tsv",
    sep="\t",
    header=True,
    index=False,
)
df_process_list[2].to_csv(path_or_buf="Table_all_details_mask_primers.tsv", sep="\t", header=True, index=False)
df_process_list[3].to_csv(path_or_buf="Table_all_details_deduplicate.tsv", sep="\t", header=True, index=False)
df_process_list[4].to_csv(path_or_buf="Table_all_details_igblast.tsv", sep="\t", header=True, index=False)

values = [
    df_process_list[0].sort_values(by=["Sample"]).iloc[:, 0].tolist(),
    df_process_list[0].sort_values(by=["Sample"]).loc[:, "start_pairs"].tolist(),
    df_process_list[0].sort_values(by=["Sample"]).loc[:, "pass_pairs"].tolist(),
    df_process_list[1].sort_values(by=["Sample"]).loc[:, "pass_pairs"].tolist(),
    df_process_list[2].sort_values(by=["Sample"]).pivot(index="Sample", columns="readtype")["pass"]["R1"].tolist(),
    df_process_list[2].sort_values(by=["Sample"]).pivot(index="Sample", columns="readtype")["pass"]["R2"].tolist(),
    df_process_list[3].sort_values(by=["Sample"]).loc[:, "unique"].tolist(),
    df_process_list[4].sort_values(by=["Sample"]).loc[:, "repres_2"].tolist(),
    df_process_list[4].sort_values(by=["Sample"]).loc[:, "pass_igblast"].tolist(),
]

final_table = dict(zip(colnames, values))
print(final_table)
df_final_table = pd.DataFrame.from_dict(final_table)
df_final_table = df_final_table.sort_values(["Sample"], ascending=[1])


# incorporating metadata
metadata = pd.read_csv("metadata.tsv", sep="\t")
metadata = metadata[metadata.columns.drop(list(metadata.filter(regex="filename")))]
logs_metadata = metadata.merge(df_final_table, left_on="sample_id", right_on="Sample")
logs_metadata = logs_metadata.drop(["Sample"], axis=1)
logs_metadata.to_csv(path_or_buf="Table_sequences_process.tsv", sep="\t", header=True, index=False)
