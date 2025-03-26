<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-airrflow_logo_dark.png">
    <img alt="nf-core/airrflow" src="docs/images/nf-core-airrflow_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/airrflow/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/airrflow/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/airrflow/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/airrflow/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/airrflow/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.2642009-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.2642009)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/airrflow)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23airrflow-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/airrflow)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)
[![AIRR compliant](https://img.shields.io/static/v1?label=AIRR-C%20sw-tools%20v1&message=compliant&color=008AFF&labelColor=000000&style=plastic)](https://docs.airr-community.org/en/stable/swtools/airr_swtools_standard.html)

## Introduction

**nf-core/airrflow** is a bioinformatics best-practice pipeline to analyze B-cell or T-cell repertoire sequencing data. The input data can be targeted amplicon bulk sequencing data of the V, D, J and C regions of the B/T-cell receptor with multiplex PCR or 5' RACE protocol, single-cell VDJ sequencing using the 10xGenomics libraries, or assembled reads (bulk or single-cell). It can also extract BCR and TCR sequences from bulk or single-cell untargeted RNAseq data. It makes use of the [Immcantation](https://immcantation.readthedocs.io) toolset as well as other AIRRseq analysis tools.

![nf-core/airrflow overview](docs/images/airrflow_workflow_overview.png)

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/airrflow/results).

## Pipeline summary

nf-core/airrflow allows the end-to-end processing of BCR and TCR bulk and single cell targeted sequencing data. Several protocols are supported, please see the [usage documentation](https://nf-co.re/airrflow/usage) for more details on the supported protocols. The pipeline has been certified as [AIRR compliant](https://docs.airr-community.org/en/stable/swtools/airr_swtools_compliant.html) by the AIRR community, which means that it is compatible with downstream analysis tools also supporting this format.

![nf-core/airrflow overview](docs/images/metro-map-airrflow.png)

1. QC and sequence assembly

- Bulk
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
- single cell
  - cellranger vdj
    - Assemble contigs
    - Annotate contigs
    - Call cells
    - Generate clonotypes

2. V(D)J annotation and filtering (bulk and single-cell)

- Assign gene segments with `IgBlast` using a germline reference (`Change-O AssignGenes`).
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
- Create germlines and define clones, repertoire analysis (`SCOPer`, `EnchantR`).
- Build lineage trees (`Dowser`, `IgphyML`, `RAxML`, `EnchantR`).

5. Repertoire analysis and reporting

- Custom repertoire analysis pipeline report (`Alakazam`).
- Aggregate QC reports (`MultiQC`).

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

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
-r <release> \
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

For common **bulk sequencing protocols** we provide pre-set profiles that specify primers, UMI length, etc for common commercially available sequencing protocols. Please check the [Supported protocol profiles](#supported-protocol-profiles) for a full list of available profiles. An example command running the NEBNext UMI protocol profile with docker containers is:

```bash
nextflow run nf-core/airrflow \
-profile nebnext_umi,docker \
--mode fastq \
--input input_samplesheet.tsv \
--outdir results
```

A typical command to run the pipeline from **single cell raw fastq files** (10X genomics) is:

```bash
nextflow run nf-core/airrflow -r dev \
-profile <docker/singularity/podman/shifter/charliecloud/conda/institute> \
--mode fastq \
--input input_samplesheet.tsv \
--library_generation_method sc_10x_genomics \
--reference_10x reference/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz \
--outdir ./results
```

A typical command to run the pipeline from **single-cell AIRR rearrangement tables or assembled bulk sequencing fasta** data is:

```bash
nextflow run nf-core/airrflow \
-r <release> \
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

nf-core/airrflow was originally written by:

- [Gisela Gabernet](https://github.com/ggabernet)
- [Susanna Marquez](https://github.com/ssnn-airr)
- [Alexander Peltzer](https://github.com/apeltzer)

We thank the following people for their extensive assistance in the development of the pipeline:

- [David Ladd](https://github.com/dladd)
- [Friederike Hanssen](https://github.com/friederikehanssen)
- [Simon Heumos](https://github.com/subwaystation)
- [Mark Polster](https://github.com/mapo9)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#airrflow` channel](https://nfcore.slack.com/channels/airrflow) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/airrflow for your analysis, please cite the article as follows:

> **nf-core/airrflow: an adaptive immune receptor repertoire analysis workflow employing the Immcantation framework**
>
> Gisela Gabernet, Susanna Marquez, Robert Bjornson, Alexander Peltzer, Hailong Meng, Edel Aron, Noah Y. Lee, Cole G. Jensen, David Ladd, Mark Polster, Friederike Hanssen, Simon Heumos, nf-core community, Gur Yaari, Markus C. Kowarik, Sven Nahnsen, Steven H. Kleinstein. (2024) PLOS Computational Biology, 20(7), e1012265. doi: [https://doi.org/10.1371/journal.pcbi.1012265](https://doi.org/10.1371/journal.pcbi.1012265). Pubmed PMID: 39058741.

The specific pipeline version using the following DOI: [10.5281/zenodo.2642009](https://doi.org/10.5281/zenodo.2642009)

Please also cite all the tools that are being used by the pipeline. An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
