#! usr/bin/bash

nextflow run nf-core/airrflow -r 5.1.1 \
-profile docker \
--mode assembled \
--genotyping \
--input genotype_samplesheet.tsv \
--outdir test_genotype_results \
-resume
