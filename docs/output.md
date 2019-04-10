# nf-core/bcellmagic: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [Fetching databases](#fetching-db) - Fetching igblast and imgt databases
* [FastQC](#fastqc) - read quality control
* [Filter sequence quality](#filter-seq-q) - filter sequences by quality
* [Mask primers](#mask-primers) - Masking primers
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline

## Fetching databases
Fetching igblast and imgt databases.

**Output directory: `results/dbs`**
If saveDBs parameter is set, then database cache will be saved in the results directory.

* `igblast_base`
    * Contains igblast database cache.
* `imgtdb_base`
    * Contains imgt database cache.

## FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

**Output directory: `results/fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

## Filter sequence quality
Filters reads that are below a quality threshold by using the tool [FilterSeq](https://presto.readthedocs.io/en/version-0.5.11/tools/FilterSeq.html) from the Presto Immcantation toolset. The default quality threshold is 20.

**Output directory: `results/filter_by_sequence_quality`**

* `command_log.txt`
    * Log of the process that will be parsed to generate a report.
* `fastq/*.fastq`
    * fastq with only reads that passed the quality filter.
* `fastq/*.tab`
    * table containing read ID and quality.
* `fastq/*.log`
    * Log of the process.

## Mask primers
Masks primers that are provided in the C-primers and V-primers input files. It uses the tool [MaskPrimers](https://presto.readthedocs.io/en/version-0.5.11/tools/MaskPrimers.html) of the Presto Immcantation toolset.

**Output directory: `results/

## Pair sequences
Pair read mates using [PairSeq](https://presto.readthedocs.io/en/version-0.5.11/tools/PairSeq.html) from the Presto Immcantation toolset.

## MultiQC
[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see http://multiqc.info
