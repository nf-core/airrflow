#!/bin/bash
nextflow run . -profile standard,docker --metadata '/Users/alexanderpeltzer/Downloads/metasheet_test.tsv' --cprimers '/Users/alexanderpeltzer/Downloads/CPrimers_IG.fasta' --vprimers '/Users/alexanderpeltzer/Downloads/VPrimers.fasta' --max_memory 8.GB --max_cpus 8 -resume -dump-channels
