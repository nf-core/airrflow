#! /usr/bin/bash

nextflow run nf-core/airrflow -r 5.1.1 \
-profile singularity \
--mode assembled \
--genotyping \
--single_clone_representative false \
--skip_clonal_analysis \
--input genotype_samplesheet.tsv \
--outdir test_genotype_results \
-c resource.config \
-resume
