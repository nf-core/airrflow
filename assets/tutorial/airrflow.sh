nextflow run nf-core/airrflow -r 4.1.0 \
-profile docker \
--mode assembled \
--input samplesheet.tsv \
--outdir results \
-w work \
--max_cpus 12 \
--max_memory 12.GB \
--skip_multiqc
