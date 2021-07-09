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

    ID	R1	R2	I1	Source	Treatment	Extraction_time	Population
    Sample1	Sample1_dn_R1.fastq.gz	Sample1_dn_R2.fastq.gz	Sample1_dn_I1.fastq.gz	Patient1	DrugTreatment	baseline	dn
    Sample2	Sample2_m_R1.fastq.gz	Sample2_m_R2.fastq.gz	Sample2_m_I1.fastq.gz	Patient1	DrugTreatment	baseline	m
    """

    sample_run_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 7
        HEADER = ["ID", "R1", "R2", "I1", "Source", "Treatment", "Extraction_time", "Population"]
        HEADER_NOI1 = ["ID", "R1", "R2", "Source", "Treatment", "Extraction_time", "Population"]
        header = [x.strip('"') for x in fin.readline().strip().split("\t")]
        if not (header[: len(HEADER)] == HEADER or header[: len(HEADER)] == HEADER_NOI1):
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            print("or  {} != {}".format(",".join(header), ",".join(HEADER_NOI1)))

            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split("\t")]

            ## Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN)


if __name__ == "__main__":
    sys.exit(main())