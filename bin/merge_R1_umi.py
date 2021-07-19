#!/usr/bin/env python3
from Bio import SeqIO
from gzip import open as gzopen
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-R1", "--R1_filename", type=str, help="Path of R1 fastq.gz file.")
parser.add_argument(
    "-I1",
    "--I1_filename",
    type=str,
    help="Path of I1 fastq.gz file containing the UMI barcode.",
)
parser.add_argument(
    "-s", "--umi_start", type=int, help="UMI start position in index file."
)
parser.add_argument("-l", "--umi_length", type=int, help="UMI length.")
parser.add_argument(
    "-o", "--outname", type=str, default="UMI_R1.fastq.gz", help="Name of output file."
)

args = parser.parse_args()

records_I1 = SeqIO.parse(gzopen(args.I1_filename, "rt"), format="fastq")
records_R1 = SeqIO.parse(gzopen(args.R1_filename, "rt"), format="fastq")
umi_start = args.umi_start
umi_end = args.umi_start + args.umi_length

i = 0

with gzopen(args.outname, "wt") as out_handle:
    for record_I1 in records_I1:
        i += 1
        record_R1 = next(records_R1)
        record_I1 = record_I1[umi_start:umi_end]
        if record_R1.id == record_I1.id:
            record_UMI_R1 = record_I1 + record_R1
            SeqIO.write(record_UMI_R1, out_handle, "fastq")
            if i % 100000 == 0:
                print(i)
        else:
            print(record_I1, record_R1)
            exit(1)
