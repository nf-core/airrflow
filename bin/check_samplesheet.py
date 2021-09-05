#!/usr/bin/env python

# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/atacseq/design.csv


import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Read nf-core/bcellmagic samplesheet file and check its contents."
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
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in):
    """
    This function checks that the samplesheet follows the following structure:

    sample_id	filename_R1	filename_R2	filename_I1	subject_id	group_name	pcr_target_locus
    Sample1	Sample1_dn_R1.fastq.gz	Sample1_dn_R2.fastq.gz	Sample1_dn_I1.fastq.gz	Patient1	baseline_dn	ig
    Sample2	Sample2_m_R1.fastq.gz	Sample2_m_R2.fastq.gz	Sample2_m_I1.fastq.gz	Patient1	baseline_m	ig
    """

    sample_run_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 6
        HEADER = [
            "sample_id",
            "filename_R1",
            "filename_R2",
            "filename_I1",
            "subject_id",
            "group_name",
            "pcr_target_locus",
        ]
        HEADER_NOI1 = [
            "sample_id",
            "filename_R1",
            "filename_R2",
            "subject_id",
            "group_name",
            "pcr_target_locus",
        ]
        header = [x.strip('"') for x in fin.readline().strip().split("\t")]
        for col in HEADER_NOI1:
            if col not in header:
                print(
                    "ERROR: Please check samplesheet header: {} ".format(
                        ",".join(header)
                    )
                )
                print("Header is missing column {}".format(col))
                print("Header must contain columns {}".format("\t".join(HEADER_NOI1)))
                sys.exit(1)

        ## Check sample entries
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


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN)


if __name__ == "__main__":
    sys.exit(main())
