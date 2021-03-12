#!/usr/bin/env python3
"""
Clean IMGT germline fasta files for IgBLAST database build
"""
import Bio
from pkg_resources import parse_version
from sys import argv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Get input and output file names
in_file = argv[1]
out_file = argv[2]

# Load sequences into memory and process them
name_set = set()
seq_list = list()
for rec in SeqIO.parse(in_file, 'fasta'):
    name = rec.description.split('|')[1]
    if name not in name_set:
        name_set.add(name)
        seq = SeqRecord(rec.seq.ungap('.').upper(), id=name, name=name, description=name)
        seq_list.append(seq)

# Overwrite file
with open(out_file, 'w') as out_handle:
    if parse_version(Bio.__version__) >= parse_version('1.71'):
        # Biopython >= v1.71
        SeqIO.write(seq_list, out_handle, format='fasta-2line')
    else:
        # Biopython < v1.71
        writer = SeqIO.FastaIO.FastaWriter(out_handle, wrap=None)
        writer.write_file(seq_list)
