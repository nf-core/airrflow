#! usr/bin/bash

nextflow run nf-core/airrflow -r 4.2.0 \
-profile docker \
--mode fastq \
--input 10x_sc_raw.tsv \
--library_generation_method sc_10x_genomics \
--reference_10x refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0 \
-c resource.config \
--clonal_threshold 0 \
--outdir sc_from_fastq_results
