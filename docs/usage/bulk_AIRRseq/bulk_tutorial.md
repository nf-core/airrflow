# nf-core/airrflow: Bulk Airr-seq tutorial

This tutorial provides a step by step introduction on how to run nf-core/airrflow on bulk Airr-seq data.

## Pre-requisites

> [!INSTALLATION]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set up Nextflow and a container engine needed to run this pipeline. At the moment, nf-core/airrflow does NOT support using conda virtual environments for dependency management, only containers are supported. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

For the purpose of running this tutorial on your local machine, we recommend a docker installation.

To install docker, follow the instructions [here](https://docs.docker.com/engine/install/). After installation on linux, don't forget to check the [post-installation steps](https://docs.docker.com/engine/install/linux-postinstall/).

## Testing the pipeline with built-in tests

Once you have set up your Nextflow and container (docker or singularity), test the airrflow pipeline with built-in test. 

```bash
nextflow run nf-core/airrflow -r 4.2.0 -profile test,docker --outdir test_results
```

## Running airrflow pipeline 



