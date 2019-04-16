# nf-core/bcellmagic

A pipeline to analyze B-cell repertoires.

[![Build Status](https://travis-ci.org/nf-core/bcellmagic.svg?branch=master)](https://travis-ci.org/nf-core/bcellmagic)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
[![DOI](https://zenodo.org/badge/160348245.svg)](https://zenodo.org/badge/latestdoi/160348245)


[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/bcellmagic.svg)](https://hub.docker.com/r/nfcore/bcellmagic)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction

The nf-core/bcellmagic pipeline is built to analyze B-cell repertoire sequencing data from targeted amplification experiments. It makes use of the [Immcantation 2.5.0](https://immcantation.readthedocs.io/en/version-2.5.0/) toolset for the analysis of B-cell repertoires.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.

## Documentation
The nf-core/bcellmagic pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Credits
nf-core/bcellmagic was originally written by Gisela Gabernet, Simon Heumos, Alexander Peltzer.
