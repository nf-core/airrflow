# nf-core/airrflow: Single-cell tutorial

This tutorial provides a step by step introduction on how to run nf-core/airrflow on single-cell BCR-seq data or single-cell TCR-seq data.

## Pre-requisites

> [!INSTALLATION]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set up Nextflow and a container engine needed to run this pipeline. At the moment, nf-core/airrflow does NOT support using conda virtual environments for dependency management, only containers are supported. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

For the purpose of running this tutorial on your local machine, we recommend a docker installation.

To install docker, follow the instructions [here](https://docs.docker.com/engine/install/). After installation on linux, don't forget to check the [post-installation steps](https://docs.docker.com/engine/install/linux-postinstall/).

## Testing the pipeline with built-in tests

Once you have set up your Nextflow and container (docker or singularity). Test the airrflow pipeline with built-in tests. 

```bash
nextflow run nf-core/airrflow -r 4.2.0 -profile test,docker --outdir test_results
```


## Running airrflow pipeline
There are two acceptable input formats for airrflow single-cell AIRRseq pipeline: AIRR rearrangement or fastq format. 

For this tutorial we will practice on both of the input formats. 

## AIRR rearrangement format
### Datasets

For this tutorial we will use subsampled PBMC single-cell BCR sequencing data from two subjects, before (d0) and after flu vaccination (d12). The dataset is available on [Zenodo](https://zenodo.org/doi/10.5281/zenodo.11373740). You don't need to downlaod the samples bacause the link is already in the samplesheet. 

### Preparing samplesheet and configuration file

To run the pipeline one needs to prepare a samplesheet and configuration file. 

Ready samplesheet for this tutorial is [here](sample_data_code/assembled_samplesheet.tsv).
Ready configuration file is [here](sample_data_code/resource.config). 
Download these two files to the directory where you will run the airrflow pipeline. 

### Running airrflow 

With all the files prepared, now you are ready to run the airrflow pipeline. 

```bash
nextflow run nf-core/airrflow -r 4.2.0 \
-profile docker \
--mode assembled \
--input assembled_samplesheet.tsv \
--outdir sc_from_assembled_results  \
-c resource.config
```
Of course you can wrap all your code in a bash file. We prepared one for you and it's available [here](sample_data_code/airrflow_sc_from_assembled.sh).
With the bash file, it's easy to run the pipeline with a single-line command. 

```bash
bash airrflow_sc_from_assembled.sh
```


## fastq format
### Datasets

### Preparing samplesheet, gene reference and configuration file

### Running airrflow


## Understanding the results

Explanation of the analysis steps of the pipeline and where to find the results. Point to output documentation for specific steps if necessary.

## Including lineage tree computation

Re-run with --lineage_tree parameter using -resume as example.

## Downstream analysis

Downstream analysis can be performed from the AIRR repertoires. Provide one example and links to the Immcantation single-cell tutorial.
