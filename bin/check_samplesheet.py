#!/usr/bin/env python
# Written by Gisela Gabernet and released under the MIT license (2020).
# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/atacseq/design.csv

import os
import sys
import errno
import argparse
import pandas as pd
import re


def parse_args(args=None):
    Description = "Read nf-core/airrflow samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("file_in", help="Input samplesheet file.")
    parser.add_argument("-a", "--assembled", help="Input samplesheet type", action="store_true", default=False)
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


def check_samplesheet(file_in, assembled):
    """
    This function checks that the samplesheet:

    - contains the compulsory fields: sample_id, filename_R1, filename_R2, subject_id, pcr_target_locus, species, single_cell
    - sample ids are unique
    - samples from the same subject come from the same species
    - pcr_target_locus is "IG" or "TR"
    - species is "human" or "mouse"
    """

    sample_run_dict = {}
    with open(file_in, "r") as fin:
        # Defining minimum columns and required columns
        min_cols = 7
        required_columns_raw = [
            "sample_id",
            "filename_R1",
            "filename_R2",
            "subject_id",
            "species",
            "pcr_target_locus",
            "single_cell",
            "sex",
            "tissue",
            "biomaterial_provider",
            "age",
        ]
        required_columns_assembled = [
            "sample_id",
            "filename",
            "subject_id",
            "species",
            "pcr_target_locus",
            "single_cell",
            "sex",
            "tissue",
            "biomaterial_provider",
            "age",
        ]
        no_whitespaces_raw = [
            "sample_id",
            "filename_R1",
            "filename_R2",
            "subject_id",
            "species",
            "pcr_target_locus",
            "tissue",
        ]
        no_whitespaces_assembled = [
            "sample_id",
            "filename",
            "subject_id",
            "species",
            "pcr_target_locus",
            "tissue",
        ]

        ## Read header
        header = [x.strip('"') for x in fin.readline().strip().split("\t")]
        ## Read tab
        tab = pd.read_csv(file_in, sep="\t", header=0)

        # Check that all required columns for assembled and raw samplesheets are there, and do not contain whitespaces
        if assembled:
            for col in required_columns_assembled:
                if col not in header:
                    print("ERROR: Please check samplesheet header: {} ".format(",".join(header)))
                    print("Header is missing column {}".format(col))
                    print("Header must contain columns {}".format("\t".join(required_columns)))
                    raise IndexError("Header must contain columns {}".format("\t".join(required_columns)))
            for col in no_whitespaces_assembled:
                values = tab[col].tolist()
                if any([re.search(r"\s+", s) for s in values]):
                    print_error(
                        "The column {} contains values with whitespaces. Please ensure that there are no tabs, spaces or any other whitespaces in these columns as well: {}".format(
                            col, no_whitespaces_assembled
                        )
                    )

        else:
            for col in required_columns_raw:
                if col not in header:
                    print("ERROR: Please check samplesheet header: {} ".format(",".join(header)))
                    print("Header is missing column {}".format(col))
                    print("Header must contain columns {}".format("\t".join(required_columns)))
                    raise IndexError("Header must contain columns {}".format("\t".join(required_columns)))
            for col in no_whitespaces_raw:
                values = tab[col].tolist()
                if any([re.search(r"\s+", s) for s in values]):
                    print_error(
                        "The column {} contains values with whitespaces. Please ensure that there are no tabs, spaces or any other whitespaces in these columns as well: {}".format(
                            col, no_whitespaces_raw
                        )
                    )

        ## Check that rows have the same fields as header, and at least the compulsory ones are provided
        for line_num, line in enumerate(fin):
            lspl = [x.strip().strip('"') for x in line.strip().split("\t")]

            ## Check valid number of columns per row
            if len(lspl) < len(header):
                print_error(
                    "Invalid number of columns in this row (should be {})!".format(len(header)),
                    "Line {}".format(line_num),
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < min_cols:
                print_error(
                    "Invalid number of populated columns (should be {})!".format(min_cols),
                    "Line",
                    line,
                )

        ## Check that sample ids are unique
        if len(tab["sample_id"]) != len(set(tab["sample_id"])):
            print_error(
                "Sample IDs are not unique! The sample IDs in the input samplesheet should be unique for each sample."
            )

        ## Check that pcr_target_locus is IG or TR
        for val in tab["pcr_target_locus"]:
            if val.upper() not in ["IG", "TR"]:
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
    check_samplesheet(args.file_in, args.assembled)


if __name__ == "__main__":
    sys.exit(main())
