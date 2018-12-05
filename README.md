# nf-core/bcellmagic

A pipeline to analyze B-cell and T-cell repertoires based on [Immcantation 2.5.0](https://immcantation.readthedocs.io/en/version-2.5.0/).

[![Build Status](https://travis-ci.org/nf-core/bcellmagic.svg?branch=master)](https://travis-ci.org/nf-core/bcellmagic)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/bcellmagic.svg)](https://hub.docker.com/r/nfcore/bcellmagic)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.


### Documentation
The nf-core/bcellmagic pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

Input needs to be a TSV file following this format in general:

```
ID	Source	Treatment	Extraction_time	Population	R1	R2	I1
QMKMK072AD	Patient 2	Terifluomid	baseline	p	sample_S8_L001_R1_001.fastq.gz	sample_S8_L001_R2_001.fastq.gz	sample_S8_L001_I1_001.fastq.gz
```
Attention, the R1/R2 and I1 naming patterns are crucial!

An example call of the pipeline could be then:

```
nextflow run ggabernet/bcellmagic -profile standard,docker --metadata metasheet_test.tsv --cprimers CPrimers_IG.fasta --vprimers VPrimers.fasta --max_memory 8.GB --max_cpus 8 -resume 
```

### Credits
nf-core/bcellmagic was originally written by Gisela Gabernet, Simon Heumos and Alexander Peltzer.
