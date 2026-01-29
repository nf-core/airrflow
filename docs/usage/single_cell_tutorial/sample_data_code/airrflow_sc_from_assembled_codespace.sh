#! /usr/bin/bash

nextflow run nf-core/airrflow -r 4.4.0 \
-profile singularity \
--mode assembled \
--input assembled_samplesheet.tsv \
--outdir sc_from_assembled_results  \
-c resource.config \
-resume
