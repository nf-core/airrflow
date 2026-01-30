# nf-core/airrflow: Single-cell tutorial

This tutorial provides a step by step introduction on how to run nf-core/airrflow on single-cell BCR-seq data or single-cell TCR-seq data.

<iframe width="560" height="315" src="https://www.youtube.com/embed/TT6eaoiZQvg?si=Psbo9cQ7KemfiBd5" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>

## Pre-requisites

You can run this tutorial using the Github Codespaces platform. Codespaces already has Nextflow and Singularity pre-installed, and it can automatically be used for every nf-core repository. To create a Codespace instance for nf-core/airrflow, first click on the button labelled `Code` at the top of [nf-core/airrflow repository](https://github.com/nf-core/airrflow).

In the dropdown menu, go to the `Codespaces` tab. Click the `...` sign, then select `+ New with options...`.

<p align="center">
<img src="https://raw.githubusercontent.com/nf-core/airrflow/dev/docs/images/Create_codespaces.png" width="400" alt="Create Codespaces with options">
</p>

After that, you’ll be directed to the configuration page. Select "4-core" for `machine type`, which will give you 4 CPUs, 16GB RAM and 32GB space.

<p align="center">
<img src="https://raw.githubusercontent.com/nf-core/airrflow/dev/docs/images/Codespaces_4core.png" width="400" alt="Chose 4-core">
</p>

If you want to know more about Codespaces, check [the Codespaces overview](https://docs.github.com/en/codespaces/about-codespaces/what-are-codespaces) or the Codespaces section in nf-core documentation [the Devcontainers overview](https://nf-co.re/docs/tutorials/devcontainer/overview).

When running this tutorial on your local machine, you'll first have to set up Nextflow and a container engine (Docker or Singularity).

> [!NOTE]
> If you want to run this tutorial on your local machine, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set up Nextflow and a container engine needed to run this pipeline. At the moment, nf-core/airrflow does NOT support using conda virtual environments for dependency management, only containers are supported. To install Docker, follow the [instructions](https://docs.docker.com/engine/install/). After installation Docker on Linux, don't forget to check the [post-installation steps](https://docs.docker.com/engine/install/linux-postinstall/).

## Testing the pipeline with built-in tests

Once you have set up Nextflow and container (Docker or Singularity) for your local machine or Codespaces environment, test nf-core/airrflow with the built-in test data.

```bash
nextflow run nf-core/airrflow -r 5.0 -profile test,docker --outdir test_results
```

Change the `docker` profile to `singularity` if you use Codespaces since Docker currently cannot be used in Codespaces. You can first set up a Singularity cache directory which will allow the reuse of Singularity container across all runs:

```bash
mkdir singularity_cache
export NXF_SINGULARITY_CACHEDIR="/workspaces/airrflow/singularity_cache"
```

Then run nf-core/airrflow with the test data:

```bash
nextflow run nf-core/airrflow -r 5.0 -profile test,singularity --outdir test_results
```

> [!NOTE]
> The '-r' flag in the command specifies which nf-core/airrflow release to run. We recommend always [checking and using the latest release] (https://nf-co.re/airrflow/releases_stats/).

> [!NOTE]
> Because Codespaces provides limited CPU and RAM resources, the test run may take 20-25 minutes. The process will be faster on systems with greater CPU and RAM capacity.

If the tests run through correctly, you should see this output in your command line:

```bash
-[nf-core/airrflow] Pipeline completed successfully-
Completed at: 25-Nov-2025 21:07:46
Duration    : 23m 56s
CPU hours   : 1.1
Succeeded   : 221
```

## Supported input formats

There are two supported input formats for nf-core/airrflow single-cell AIRR-seq pipeline: assembled sequences in AIRR rearrangement format or raw reads in fastq format sequenced in the 10x Genomics platform.

The [AIRR rearrangement format](https://docs.airr-community.org/en/latest/datarep/rearrangements.html) is a standard format to store BCR and TCR sequence data with relevant metadata fields. This format is supported as input and output by multiple tools specific for analyzing AIRR-seq data. For example, when analyzing single-cell AIRR sequencing data with CellRanger versions >= 4.0, an AIRR rearrangement file will be provided as output, and this is the recommended input for running nf-core/airrflow. Note that it is also possible to start running the pipeline directly from raw sequencing reads, and in this case CellRanger will be run when launching nf-core/airrflow.

The AIRR rearrangement format is also the default one when analyzing publicly available data from specialized AIRR-seq databases such as the AIRR Data Commons through the [iReceptor gateway](https://gateway.ireceptor.org/login).

In this tutorial we will showcase how to run nf-core/airrflow with both of the input formats.

![nf-core/airrflow overview](https://raw.githubusercontent.com/nf-core/airrflow/master/docs/images/airrflow_workflow_overview.png)

## Starting from AIRR rearrangement format

### Datasets

For this tutorial we will use sub-sampled PBMC single-cell BCR sequencing data from two subjects, before (d0) and after flu vaccination (d12).
The dataset is publicly available on [Zenodo](https://zenodo.org/doi/10.5281/zenodo.11373740).
You don't need to download the dataset because the links to the samples are already provided in the samplesheet and Nextflow will get the data from the links automatically when running the pipeline.

### Preparing the samplesheet and configuration file

To run the pipeline, a tab-separated samplesheet that provides the path to the AIRR rearrangement files must be prepared.
The samplesheet collects experimental details that are important for the data analysis.

Details on the required columns of a samplesheet are available [here](https://nf-co.re/airrflow/usage#assembled-input-samplesheet-bulk-or-single-cell-sequencing).

The resource configuration file sets the compute infrastructure maximum available number of CPUs, RAM memory and running time. This will ensure that no pipeline process requests more resources than available in the compute infrastructure where the pipeline is running. The resource config should be provided with the `-c` option. In this example we set the maximum RAM memory to 15GB, we restrict the pipeline to use 4 CPUs and to run for a maximum of 24 hours.

```json title="resource.config"
process {
    resourceLimits = [ memory: 15.GB, time: 24.h, cpus: 4 ]
}
```

We prepared the [samplesheet](https://github.com/nf-core/airrflow/blob/dev/docs/usage/single_cell_tutorial/sample_data_code/assembled_samplesheet.tsv) and the [configuration file](https://github.com/nf-core/airrflow/tree/dev/docs/usage/single_cell_tutorial/sample_data_code/resource.config) for this tutorial. Download both files to the directory where you intend to run nf-core/airrflow.

> [!TIP]
> Before setting memory and CPUs in the configuration file, we recommend verifying the available memory and CPUs on your system. Otherwise, exceeding the system's capacity may result in an error indicating that you requested more CPUs than available or run out of memory. You can also remove the "time" parameter from the configuration file to allow for unlimited runtime for large-size dataset.

> [!NOTE]
> When running nf-core/airrflow with your own data, provide the full path to your input files under the filename column.

### Running airrflow

With all the files ready, you can start the pipeline with the following command if you run it locally.

```bash
nextflow run nf-core/airrflow -r 5.0 \
-profile docker \
--mode assembled \
--input assembled_samplesheet.tsv \
--outdir sc_from_assembled_results  \
-c resource.config \
-resume
```

Of course you can wrap all your code in a [bash file](https://github.com/nf-core/airrflow/tree/dev/docs/usage/single_cell_tutorial/sample_data_code/airrflow_sc_from_assembled.sh). With the bash file, it's easy to run the pipeline with a single-line command.

```bash
bash airrflow_sc_from_assembled.sh
```

If you are running the pipeline on Codespace, remember to replace `docker` profile with `singularity`.

Or run this [bash file](https://github.com/nf-core/airrflow/tree/dev/docs/usage/single_cell_tutorial/sample_data_code/airrflow_sc_from_assembled_codespace.sh) within the folder where it locates.

```bash
bash airrflow_sc_from_assembled_codespace.sh
```

> [!TIP]
> When launching a Nextflow pipeline with the `-resume` option, any processes that have already been run with the exact same code, settings and inputs will be cached and the pipeline will resume from the last step that changed or failed with an error. The benefit of using "resume" is to avoid duplicating previous work and save time when re-running a pipeline.
> We include "resume" in our Nextflow command as a precaution in case anything goes wrong during execution. After fixing the issue, you can relaunch the pipeline with the same command, it will resume running from the point of failure, significantly reducing runtime and resource usage.

After launching the pipeline the following will be printed to the console output, followed by some the default parameters used by the pipeline and execution log of airrflow processes:

```bash
 N E X T F L O W   ~  version 24.10.5

WARN: It appears you have never run this project before -- Option `-resume` is ignored
Launching `https://github.com/nf-core/airrflow` [boring_heyrovsky] DSL2 - revision: d91dd840f4 [5.0]


------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/airrflow 5.0
------------------------------------------------------

```

Once the pipeline has finished successfully, the following message will appear:

```bash
-[nf-core/airrflow] Pipeline completed successfully-
Completed at: 17-Nov-2025 19:20:17
Duration    : 7m 10s
CPU hours   : 0.3
Succeeded   : 44
```

## Starting from raw reads in fastq format

### Datasets

For this tutorial we will use subsampled blood single-cell TCR sequencing data of one subject generated from the 10x Genomic platform. The links to the fastq files are in the samplesheet.

### Preparing samplesheet, gene reference and configuration file

To run nf-core/airrflow on single cell TCR or BCR sequencing data from fastq files, we need to prepare samplesheet, pre-built 10x genomics V(D)J references and configuration file in advance. Details on the required columns for this samplesheet are available [here](https://nf-co.re/airrflow/usage#fastq-input-samplesheet-single-cell-sequencing).

> [!WARNING]
> The fastq file names must follow the 10X Genomics file naming convention or the cellranger process will fail.

We prepared the [samplesheet](https://github.com/nf-core/airrflow/blob/dev/docs/usage/single_cell_tutorial/sample_data_code/10x_sc_raw.tsv) and the [configuration file](https://github.com/nf-core/airrflow/blob/dev/docs/usage/single_cell_tutorial/sample_data_code/resource.config) for this tutorial. Download these two files to the directory where you intend to run nf-core/airrflow.

> [!TIP]
> Before setting memory and CPUs in the configuration file, we recommend verifying the available memory and CPUs on your system. Otherwise, exceeding the system's capacity may result in an error indicating that you requested more CPUs than available or run out of memory.

Pre-built 10x genomics V(D)J references can be accessed at the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/downloads). Both human and mouse V(D)J references are available. Download the reference that corresponds to the species of your dataset.

### Running airrflow

To run the pipeline locally, use the following command to launch nf-core/airrflow for the dataset in this tutorial:

```bash
nextflow run nf-core/airrflow -r 5.0 \
-profile docker \
--mode fastq \
--input 10x_sc_raw.tsv \
--library_generation_method sc_10x_genomics \
--reference_10x refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0 \
-c resource.config \
--clonal_threshold 0 \
--outdir sc_from_fastq_results \
-resume
```

In this tutorial, since the samples are TCRs, which do not have somatic hypermutation, clones are defined strictly by identical junction regions. For this reason, we set the `--clonal_threshold` parameter to 0. For more details on important considerations when performing clonal analysis check [FAQ](./FAQ.md).

Of course you can wrap all your code in a bash file. We prepared one for you and it's available [here](https://github.com/nf-core/airrflow/tree/dev/docs/usage/single_cell_tutorial/sample_data_code/airrflow_sc_from_fastq.sh).
With the bash file, it's easy to run the pipeline with a single-line command.

```bash
bash airrflow_sc_from_fastq.sh
```

> [!NOTE]
> Due to the limited RAM and storage space, the single cell raw reads example in this tutorial cannot be run on Codespace.

After launching the pipeline the following will be printed to the console output, followed by some the default parameters used by the pipeline and execution log of airrflow processes:

```bash
 N E X T F L O W   ~  version 24.10.5

Launching `https://github.com/nf-core/airrflow` [gloomy_monod] DSL2 - revision: d91dd840f4 [5.0]


------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/airrflow 5.0
------------------------------------------------------
```

Once the pipeline has finished successfully, the following message will appear:

```bash
-[nf-core/airrflow] Pipeline completed successfully-
Completed at: 11-Mar-2025 13:18:13
Duration    : 2m 46s
CPU hours   : 0.3 (0.1% cached)
Succeeded   : 17
Cached      : 2
```

## Understanding the results

After running the pipeline, several sub-folders are available under the results folder.

```bash
Airrflow_report.html
- cellranger
- vdj_annotation
- qc_filtering
- clonal_analysis
- repertoire_comparison
- multiqc
- report_file_size
- pipeline_info
```

The summary report, named `Airrflow_report.html`, provides an overview of the analysis results, such as an overview of the number of sequences per sample in each of the pipeline steps, the V(D)J gene assignment and QC, and V gene family usage. Additionally, it contains links to detailed reports for other specific analysis steps.

The analysis steps and their corresponding folders, where the results are stored, are briefly listed below. Detailed documentation on the pipeline output can be found on the [Output documentation page](https://nf-co.re/airrflow/docs/output/).

1. QC and sequence assembly (if starting from fastq files).
   - In this first step, Cell Ranger's VDJ algorithm is employed to assemble contigs, annotate contigs, and call cells. The results are stored in the 'cellranger' folder.

2. V(D)J annotation and filtering.
   - In this step, V(D)J gene segments are inferred using the provided germline reference and [`IgBLAST`](https://www.ncbi.nlm.nih.gov/igblast/). Alignments are annotated in AIRR format. Non-productive sequences and sequences with low alignment quality are filtered out unless otherwise specified. The intermediate results are stored under the folder named 'vdj_annotation'.

3. QC filtering.
   - In this step, cells without heavy chains or with multiple heavy chains are removed. Sequences in different samples that share both the same cell_id and nucleotide sequence are filtered out. The result are stored in the 'qc-filtering' folder.

4. Clonal analysis.
   - Results of the clonal threshold determination using `SHazaM` should be inspected in the html report under the 'clonal_analysis/find_threshold' folder. If the automatic threshold is unsatisfactory, you can set the threshold manually and re-run the pipeline.
     (Tip: use -resume whenever running the Nextflow pipeline to avoid duplicating previous work).
   - Clonal inference is performed with `SCOPer`. Clonal inference results as well as clonal abundance and diversity plots can be inspected in the html report in the folder 'clonal_analysis/clonal_assignment'. For BCR sequencing data, mutation frequency is also computed using `SHazaM` at this step and plotted in the report. The `repertoires` subfolder contains the AIRR formatted files with the clonal assignments in a new column `clone_id` and mutation frequency in the column `mu_freq`. The `tables` subfolder contains the tabulated abundance and diversity computation as well as a table with the number of clones and their size. The `ggplots` subfolder contains the abundance and diversity plots as an `RData` object for loading and customization in R.
   - If lineage trees were computed using `Dowser`, a folder under 'clonal_analysis/dowser_lineages' will be present. The trees can be inspected in the html report and saved as PDF. Additionally, an `RDS` object with the formatted trees can also be loaded in R for customizing the lineage tree plots with Dowser.

5. Repertoire analysis
   - Comparison of several repertoire characteristics, such as V gene usage, across subjects, time points or cell populations. All associated plots and tables are available under the folder `repertoire_comparison`. The plots are also included in the `Airrflow_report.html` file. This report is generated from an R markdown `Rmd` file. It is possible to customize this to meet the user's needs by editing the report and then providing the edited Rmd file with the `--report_rmd` parameter. Check the remaining [Report parameters](https://nf-co.re/airrflow/parameters/#report-options) for further customizing the report.

6. Other reporting.
   Additional reports are also generated, including:
   - MultiQC report: summarizes QC metrics across all samples.
   - Pipeline_info report: various reports relevant to the running and execution of the pipeline.
   - Report_file_size report: Summary of the number of sequences left after each of the most important pipeline steps.

## Find out more

To continue learning about how to use nf-core/airrflow please check out the following documentation:

- [nf-core/airrflow usage documentation](https://nf-co.re/airrflow/docs/usage)
- [nf-core/airrflow parameters documentation](https://nf-co.re/airrflow/parameters)
- [FAQ page](./FAQ.md)

The nf-core troubleshooting documentation will also help you troubleshoot your Nextflow errors

- [nf-core troubleshooting](https://nf-co.re/docs/usage/troubleshooting/overview)
