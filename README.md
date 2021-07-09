# ![nf-core/bcellmagic](docs/images/nf-core-bcellmagic_logo.png)

**B cell repertoire analysis pipeline with the Immcantation framework.**.

[![GitHub Actions CI Status](https://github.com/nf-core/bcellmagic/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/bcellmagic/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/bcellmagic/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/bcellmagic/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A521.04.1-brightgreen.svg)](https://www.nextflow.io/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3607408.svg)](https://doi.org/10.5281/zenodo.3607408)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/bcellmagic.svg)](https://hub.docker.com/r/nfcore/bcellmagic)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23bcellmagic-4A154B?logo=slack)](https://nfcore.slack.com/channels/bcellmagic)

## Introduction

The nf-core/bcellmagic pipeline is built to analyze B-cell or T-cell bulk repertoire sequencing data. It makes use of the [Immcantation](https://immcantation.readthedocs.io) toolset and requires targeted amplicon sequencing data of the V, D, J and C regions of the B/T-cell receptor with multiplex PCR or 5' RACE protocol.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation) (`>=21.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/bcellmagic -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    ```bash
    nextflow run nf-core/bcellmagic -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input "samplesheet.tsv" --protocol "pcr_umi" --cprimers "CPrimers.fasta" --vprimers "VPrimers.fasta" --umi_length 12 --loci "ig"
    ```

See [usage docs](https://nf-co.re/bcellmagic/usage) for all of the available options when running the pipeline.

## Pipeline Summary

By default, the pipeline currently performs the following steps:

* Raw read quality control (`FastQC`)
* Preprocessing (`pRESTO`)
  * Filtering sequences by sequencing quality.
  * Masking amplicon primers.
  * Pairing read mates.
  * Cluster sequences according to similarity, it helps identify if the UMI diversity was not high enough.
  * Building consensus of sequences with the same UMI barcode.
  * Re-pairing read mates.
  * Assembling R1 and R2 read mates.
  * Removing and annotating read duplicates with different UMI barcodes.
  * Filtering out sequences that do not have at least 2 duplicates.
* Assigning gene segment alleles from teh IgBlast database (`Change-O`).
* Determining the BCR / TCR genotype of the sample and finding the threshold for clone definition (`TIgGER`, `SHazaM`).
* Clonal assignment: defining clonal lineages of the B-cell / T-cell populations (`Change-O`).
* Reconstructing gene calls of germline sequences (`Change-O`).
* Generating clonal trees (`Alakazam`).
* Clonal analysis (`Alakazam`).
* Repertoire comparison: calculation of clonal diversity and abundance (`Alakazam`).
* Aggregating QC reports (`MultiQC`).

## Documentation

The nf-core/bcellmagic pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/bcellmagic/usage) and [output](https://nf-co.re/bcellmagic/output).

## Credits

nf-core/bcellmagic was originally written by Gisela Gabernet, Simon Heumos, Alexander Peltzer.

<!-- We thank the following people for their extensive assistance in the development
of this pipeline: -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#bcellmagic` channel](https://nfcore.slack.com/channels/bcellmagic) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use  nf-core/bcellmagic for your analysis, please cite it using the following DOI: [10.5281/zenodo.3607408](https://doi.org/10.5281/zenodo.3607408)

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

In addition, references of tools and data used in this pipeline are as follows:

* **pRESTO** Vander Heiden, J. A., Yaari, G., Uduman, M., Stern, J. N. H., O’Connor, K. C., Hafler, D. A., … Kleinstein, S. H. (2014). pRESTO: a toolkit for processing high-throughput sequencing raw reads of lymphocyte receptor repertoires. Bioinformatics, 30(13), 1930–1932. [https://doi.org/10.1093/bioinformatics/btu138](https://doi.org/10.1093/bioinformatics/btu138).
* **SHazaM, Change-O** Gupta, N. T., Vander Heiden, J. A., Uduman, M., Gadala-Maria, D., Yaari, G., & Kleinstein, S. H. (2015). Change-O: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data: Table 1. Bioinformatics, 31(20), 3356–3358. [https://doi.org/10.1093/bioinformatics/btv359](https://doi.org/10.1093/bioinformatics/btv359).
* **Alakazam** Stern, J. N. H., Yaari, G., Vander Heiden, J. A., Church, G., Donahue, W. F., Hintzen, R. Q., … O’Connor, K. C. (2014). B cells populating the multiple sclerosis brain mature in the draining cervical lymph nodes. Science Translational Medicine, 6(248). [https://doi.org/10.1126/scitranslmed.3008879](https://doi.org/10.1126/scitranslmed.3008879).
* **TIgGER** Gadala-maria, D., Yaari, G., Uduman, M., & Kleinstein, S. H. (2015). Automated analysis of high-throughput B-cell sequencing data reveals a high frequency of novel immunoglobulin V gene segment alleles. Proceedings of the National Academy of Sciences, 112(8), 1–9. [https://doi.org/10.1073/pnas.1417683112](https://doi.org/10.1073/pnas.1417683112).
* **FastQC** Download: [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* **MultiQC** Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. [https://doi.org/10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354). Download: [https://multiqc.info/](https://multiqc.info/).
