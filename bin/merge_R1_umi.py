#!/usr/bin/env python2
from Bio import SeqIO
import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-R1", "--R1_filename", type=str, help="Path of R1 fastq.gz file.")
parser.add_argument("-I1", "--I1_filename", type=str, help="Path of I1 fastq.gz file containing the UMI barcode.")
parser.add_argument("-o", "--outname", type=str, default="UMI_R1.fastq.gz", help="Name of output file.")

args = parser.parse_args()

records_I1 = SeqIO.parse(gzip.open(args.I1_filename, "rt"), format="fastq")
records_R1 = SeqIO.parse(gzip.open(args.R1_filename, "rt"), format="fastq")

i = 0

with gzip.open(args.outname, "wt") as out_handle:
    for record_I1 in records_I1:
        i += 1
        record_R1 = next(records_R1)
        record_I1 = record_I1[6:14]
        if (record_R1.id == record_I1.id):
            record_UMI_R1 = record_I1 + record_R1
            SeqIO.write(record_UMI_R1, out_handle, "fastq")
            if i % 100000 == 0:
                print(i)
        else:
            print(record_I1, record_R1)
            exit(1)