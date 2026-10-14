#! /usr/bin/bash

nextflow run nf-core/airrflow -r 5.1.2 \
-profile docker \
--mode assembled \
--input assembled_samplesheet.tsv \
--outdir sc_from_assembled_results  \
-c resource.config \
-resume
