nextflow run nf-core/airrflow -r 4.2.0 \
-profile docker \
--mode assembled \
--input samplesheet.tsv \
--outdir results \
-w work \
--skip_multiqc
