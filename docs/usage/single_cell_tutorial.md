# nf-core/airrflow: Single-cell tutorial

This tutorial provides a step by step introduction on how to run nf-core/airrflow on single-cell BCR-seq data.

## Pre-requisites

> [!INSTALLATION]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set up Nextflow and a container engine needed to run this pipeline. At the moment, nf-core/airrflow does NOT support using conda virtual environments for dependency management, only containers are supported. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

For the purpose of running this tutorial on your local machine, we recommend a docker installation.

To install docker, follow the instructions...

## Testing the pipeline with built-in tests

Command to test the pipeline with built-in tests.

```
nextflow run nf-core/airrflow -r 4.2.0 -profile test,docker
```

## Datasets

For this tutorial we will employ the following dataset hosted in Zenodo.

Explain the datasets content.

## Preparing the samplesheet

To run the pipeline one needs to prepare a samplesheet...

Ready samplesheet for this tutorial is [here](../../assets/single_cell_tutorial/samplesheet.tsv).

## Running airrflow

Example command to run airrflow:

```bash

nextflow run ...

```

## Understanding the results

Explanation of the analysis steps of the pipeline and where to find the results. Point to output documentation for specific steps if necessary.

## Including lineage tree computation

Re-run with --lineage_tree parameter using -resume as example. Add

## Downstream analysis

Downstream analysis can be performed from the AIRR repertoires. Provide one example and links to the Immcantation single-cell tutorial.
