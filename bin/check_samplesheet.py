#!/usr/bin/env python

# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/atacseq/design.csv

import os
import sys
import errno
import argparse
import pandas as pd

def parse_args(args=None):
    Description = "Read nf-core/airrflow samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check input samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check input samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in):
    """
    This function checks that the samplesheet:

    - contains the compulsory fields: sample_id, filename_R1, filename_R2, subject_id, pcr_target_locus, species
    - sample ids are unique
    - samples from the same subject come from the same species
    - pcr_target_locus is "IG" or "TR"
    - species is "human" or "mouse"
    """

    sample_run_dict = {}
    with open(file_in, "r") as fin:

        ## Check that required columns are present
        MIN_COLS = 6
        REQUIRED_COLUMNS = [
            "sample_id",
            "filename_R1",
            "filename_R2",
            "subject_id",
            "species",
            "pcr_target_locus",
        ]
        header = [x.strip('"') for x in fin.readline().strip().split("\t")]
        for col in REQUIRED_COLUMNS:
            if col not in header:
                print(
                    "ERROR: Please check samplesheet header: {} ".format(
                        ",".join(header)
                    )
                )
                print("Header is missing column {}".format(col))
                print(
                    "Header must contain columns {}".format("\t".join(REQUIRED_COLUMNS))
                )
                sys.exit(1)

        ## Check that rows have the same fields as header, and at least the compulsory ones are provided
        for line_num, line in enumerate(fin):
            lspl = [x.strip().strip('"') for x in line.strip().split("\t")]

            ## Check valid number of columns per row
            if len(lspl) < len(header):
                print_error(
                    "Invalid number of columns in this row (should be {})!".format(
                        len(header)
                    ),
                    "Line {}".format(line_num),
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (should be {})!".format(
                        MIN_COLS
                    ),
                    "Line",
                    line,
                )

        ## Check that sample ids are unique
        tab = pd.read_csv(file_in, sep="\t", header=0)
        if len(tab["sample_id"]) != len(set(tab["sample_id"])):
            print_error(
                "Sample IDs are not unique! The sample IDs in the input samplesheet should be unique for each sample."
            )

        ## Check that pcr_target_locus is IG or TR
        for val in tab["pcr_target_locus"]:
            if val not in ["IG", "TR"]:
                print_error("pcr_target_locus must be one of: IG, TR.")

        ## Check that species is human or mouse
        for val in tab["species"]:
            if val not in ["human", "mouse"]:
                print_error(
                    "species must be one of: human, mouse. Currently, only human or mouse reference files are supported."
                )

        ## Check that samples from the same subject are the same species
        subj_group = tab.groupby("subject_id")
        for i in range(len(subj_group)):
            if len(tab.groupby("subject_id")["species"].unique()[i]) > 1:
                print_error(
                    "The same subject_id cannot belong to different species! Check input file columns 'subject_id' and 'species'."
                )


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN)


if __name__ == "__main__":
    sys.exit(main())
