# ![nf-core/airrflow](docs/images/nf-core-airrflow_logo_light.png#gh-light-mode-only) ![nf-core/airrflow](docs/images/nf-core-airrflow_logo_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/nf-core/airrflow/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/airrflow/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/airrflow/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/airrflow/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/airrflow/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.2642009-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.2642009)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/airrflow)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23airrflow-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/airrflow)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/airrflow** is a bioinformatics best-practice pipeline to analyze B-cell or T-cell repertoire sequencing data. It makes use of the [Immcantation](https://immcantation.readthedocs.io) toolset. The input data can be targeted amplicon bulk sequencing data of the V, D, J and C regions of the B/T-cell receptor with multiplex PCR or 5' RACE protocol, or assembled reads (bulk or single cell).

![nf-core/airrflow overview](docs/images/airrflow_workflow_overview.png)

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/airrflow/results).

## Pipeline summary

nf-core/airrflow allows the end-to-end processing of BCR and TCR bulk and single cell targeted sequencing data. Several protocols are supported, please see the [usage documentation](https://nf-co.re/airrflow/usage) for more details on the supported protocols.

![nf-core/airrflow overview](docs/images/metro-map-airrflow.png)

1. QC and sequence assembly (bulk only)

- Raw read quality control, adapter trimming and clipping (`Fastp`).
- Filter sequences by base quality (`pRESTO FilterSeq`).
- Mask amplicon primers (`pRESTO MaskPrimers`).
- Pair read mates (`pRESTO PairSeq`).
- For UMI-based sequencing:
  - Cluster sequences according to similarity (optional for insufficient UMI diversity) (`pRESTO ClusterSets`).
  - Build consensus of sequences with the same UMI barcode (`pRESTO BuildConsensus`).
- Assemble R1 and R2 read mates (`pRESTO AssemblePairs`).
- Remove and annotate read duplicates (`pRESTO CollapseSeq`).
- Filter out sequences that do not have at least 2 duplicates (`pRESTO SplitSeq`).

2. V(D)J annotation and filtering (bulk and single-cell)

- Assign gene segments with `IgBlast` using the IMGT database (`Change-O AssignGenes`).
- Annotate alignments in AIRR format (`Change-O MakeDB`)
- Filter by alignment quality (locus matching v_call chain, min 200 informative positions, max 10% N nucleotides)
- Filter productive sequences (`Change-O ParseDB split`)
- Filter junction length multiple of 3
- Annotate metadata (`EnchantR`)

3. QC filtering (bulk and single-cell)

- Bulk sequencing filtering:
  - Remove chimeric sequences (optional) (`SHazaM`, `EnchantR`)
  - Detect cross-contamination (optional) (`EnchantR`)
  - Collapse duplicates (`Alakazam`, `EnchantR`)
- Single-cell QC filtering (`EnchantR`)
  - Remove cells without heavy chains.
  - Remove cells with multiple heavy chains.
  - Remove sequences in different samples that share the same `cell_id` and nucleotide sequence.
  - Modify `cell_id`s to ensure they are unique in the project.

4. Clonal analysis (bulk and single-cell)

- Find threshold for clone definition (`SHazaM`, `EnchantR`).
- Create germlines and define clones, repertoire analysis (`Change-O`, `EnchantR`).
- Build lineage trees (`SCOPer`, `IgphyML`, `EnchantR`).

5. Repertoire analysis and reporting

- Custom repertoire analysis pipeline report (`Alakazam`).
- Aggregate QC reports (`MultiQC`).

## Usage

:::note
If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
with `-profile test` before running the workflow on actual data.
:::

First, ensure that the pipeline tests run on your infrastructure:

```bash
nextflow run nf-core/airrflow -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute> --outdir <OUTDIR>
```

To run nf-core/airrflow with your data, prepare a tab-separated samplesheet with your input data. Depending on the input data type (bulk or single-cell, raw reads or assembled reads) the input samplesheet will vary. Please follow the [documentation on samplesheets](https://nf-co.re/airrflow/usage#input-samplesheet) for more details. An example samplesheet for running the pipeline on bulk BCR / TCR sequencing data in fastq format looks as follows:

| sample_id | filename_R1                     | filename_R2                     | filename_I1                     | subject_id | species | pcr_target_locus | tissue | sex    | age | biomaterial_provider | single_cell | intervention   | collection_time_point_relative | cell_subset  |
| --------- | ------------------------------- | ------------------------------- | ------------------------------- | ---------- | ------- | ---------------- | ------ | ------ | --- | -------------------- | ----------- | -------------- | ------------------------------ | ------------ |
| sample01  | sample1_S8_L001_R1_001.fastq.gz | sample1_S8_L001_R2_001.fastq.gz | sample1_S8_L001_I1_001.fastq.gz | Subject02  | human   | IG               | blood  | NA     | 53  | sequencing_facility  | FALSE       | Drug_treatment | Baseline                       | plasmablasts |
| sample02  | sample2_S8_L001_R1_001.fastq.gz | sample2_S8_L001_R2_001.fastq.gz | sample2_S8_L001_I1_001.fastq.gz | Subject02  | human   | TR               | blood  | female | 78  | sequencing_facility  | FALSE       | Drug_treatment | Baseline                       | plasmablasts |

Each row represents a sample with fastq files (paired-end).

A typical command to run the pipeline from **bulk raw fastq files** is:

```bash
nextflow run nf-core/airrflow \
-profile <docker/singularity/podman/shifter/charliecloud/conda/institute> \
--mode fastq \
--input input_samplesheet.tsv \
--library_generation_method specific_pcr_umi \
--cprimers CPrimers.fasta \
--vprimers VPrimers.fasta \
--umi_length 12 \
--umi_position R1 \
--outdir ./results
```

A typical command to run the pipeline from **single-cell AIRR rearrangement tables or assembled bulk sequencing fasta** data is:

```bash
nextflow run nf-core/airrflow \
-profile <docker/singularity/podman/shifter/charliecloud/conda/institute> \
--input input_samplesheet.tsv \
--mode assembled \
--outdir results
```

See the [usage documentation](https://nf-co.re/airrflow/usage) and the [parameter documentation](https://nf-co.re/airrflow/parameters) for more details on how to use the pipeline and all the available parameters.

:::warning
Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).
:::

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/airrflow/usage) and the [parameter documentation](https://nf-co.re/airrflow/parameters).

## Pipeline output

To see the the results of a test run with a full size dataset refer to the [results](https://nf-co.re/airrflow/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/airrflow/output).

## Credits

nf-core/airrflow was written by [Gisela Gabernet](https://github.com/ggabernet), [Susanna Marquez](https://github.com/ssnn-airr), [Alexander Peltzer](@apeltzer) and [Simon Heumos](@subwaystation).

Further contributors to the pipeline are:

- [@dladd](https://github.com/dladd)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#airrflow` channel](https://nfcore.slack.com/channels/airrflow) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/airrflow for your analysis, please cite it using the following DOI: [10.5281/zenodo.2642009](https://doi.org/10.5281/zenodo.2642009)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
