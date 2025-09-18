# nf-core/airrflow: Bulk AIRR-seq tutorial

This tutorial provides a step by step introduction on how to run nf-core/airrflow on bulk AIRR-seq data.

## Pre-requisites

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set up Nextflow and a container engine needed to run this pipeline. At the moment, nf-core/airrflow does NOT support using conda virtual environments for dependency management, only containers are supported. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) before running the workflow on actual data.

For the purpose of running this tutorial on your local machine, we recommend a Docker installation. To install Docker, follow the instructions [here](https://docs.docker.com/engine/install/). After installation Docker on Linux, don't forget to check the [post-installation steps](https://docs.docker.com/engine/install/Linux-postinstall/).

Alternatively, you can run this tutorial using the Gitpod platform which has Nextflow, nf-core and Docker pre-installed. There are three ways to open the Gitpod platform. Please watch this [video](https://www.youtube.com/watch?v=ij1msCffQZA&list=PL3TSF5whlprXVp-7Br2oKwQgU4bji1S7H&index=2) to set it up. If you want to know more about Gitpod, check [the Gitpod overview](https://nf-co.re/docs/tutorials/gitpod/overview).

## Testing the pipeline with built-in tests

Once you have set up your Nextflow and container (Docker or Singularity), test nf-core/airrflow with the built-in test data.

```bash
nextflow run nf-core/airrflow -r 4.3.0 -profile test,docker --outdir test_results
```

If the tests run through correctly, you should see this output in your command line:

```bash
-[nf-core/airrflow] Pipeline completed successfully-
Completed at: 11-Mar-2025 11:30:35
Duration    : 5m 50s
CPU hours   : 0.6
Succeeded   : 221
```

## Datasets

In this tutorial, we will use nf-core/airrflow to analyze bulk BCR sequencing data from two subjects with multiple sclerosis publicly available on Sequence Read Archive (PRJNA248475) from [Stern et al](https://pubmed.ncbi.nlm.nih.gov/25100741/). The first subject has 3 samples available from a cervical lymph node and 1 sample from a brain lesion while the second subject has 3 samples from the lymph node and 3 samples from a brain lesion. You don't need to download the samples because the links to the samples are already provided in the samplesheet.

## Preparing the samplesheet and configuration file

To run the pipeline on bulk BCR/TCR sequencing data, several files must be prepared in advance:

- A tab-separated samplesheet containing the information of each sample. Details on the required columns of a samplesheet are available [here](https://nf-co.re/airrflow/usage#input-samplesheet).
- A configuration file specifying the system's maximum available RAM memory, CPUs and running time. This will ensure that no pipeline process requests more resources than available in the compute infrastructure where the pipeline is running. The resource configuration file is provided with the `-c` option. In this example we set the maximum RAM memory to 20GB, we restrict the pipeline to use 4 CPUs and to run for a maximum of 24 hours. Depending on the size of your dataset, it might be required to extend the running time. You can also remove the "time" parameter from the configuration file to allow for unlimited runtime.

```json title="resource.config"
process {
   resourceLimits = [cpus: 4, memory: 20.GB, time: 24.h]
}
```

> [!TIP]
> Before setting memory and cpus in the configuration file, we recommend verifying the available memory and cpus on your system. Otherwise, exceeding the system's capacity may result in an error indicating that you requested more cpus than available or run out of memory.

> [!NOTE]
> When running nf-core/airrflow with your own data, provide the full path to your input files under the filename column.

A prepared samplesheet for this tutorial can be found [here](https://github.com/nf-core/airrflow/tree/dev/docs/usage/bulk_tutorial/bulk_sample_code/metadata_pcr_umi_airr_300.tsv), and the configuration file is available [here](https://github.com/nf-core/airrflow/tree/dev/docs/usage/bulk_tutorial/bulk_sample_code/resource.config).
Download both files to the directory where you intend to run nf-core/airrflow.

## Choosing the right protocol profile

Bulk BCR and TCR targeted sequencing can be performed with a wide variety of protocols, using different library preparation methods. Different protocols usually use different amplification primers, UMI barcode lengths and position, which require different parameter settings to run nf-core/airrflow. To make it easier to run the pipeline on commonly used commercially available kits, we provide parameter presets as profiles. A full [list of protocol profiles](https://nf-co.re/airrflow/docs/usage/#supported-protocol-profiles) is available on the usage documentation page.

You can provide a protocol profile with the `-profile` parameter, followed by other profiles, such as the container engine profile in a comma separated fashion. You will then usually only need to provide the input samplesheet, resource config file and output directory path. However, if you want to override any option or add additional parameters, you can provide them to the airrflow launching command as any parameters in the launch command will override the parameters in the profile.

```bash
nextflow run nf-core/airrflow -r 4.3.0 \
-profile <protocol-profile-name>,docker \
--input samplesheet.tsv \
-c resource.config \
--outdir bulk_fastq_results \
-resume
```

> [!TIP]
> We're always looking forward to expanding the set of protocol profiles readily available for other users. Feel free to open an issue and create a pull request to add a new profile that you want to share with other users or ask in the nf-core `#airrflow` [slack channel](https://nf-co.re/join) if you have any questions in doing so.

## Analyzing a dataset with a custom library preparation method

If your dataset was generated using a custom library preparation method, you can manually set the relevant parameters according to your protocol design, similar to the approach we used for the samples in this tutorial. For more examples on how to set the parameters for custom protocols check the [usage documentation](https://nf-co.re/airrflow/docs/usage/#supported-bulk-library-generation-methods-protocols) page.

The BCRseq dataset used in this tutorial was obtained with a multiplexed PCR protocol using custom C-region and V-region primers. We stored the sequences for the V-region primers as well as the C-region primers in AWS S3, and the links are provided in the Nextflow command which will be fetched by Nextflow automatically when executing the command. You can also provide the full path to the custom primers fasta files.

The command to launch nf-core/airrflow for the dataset in this tutorial is the following:

```bash
nextflow run nf-core/airrflow -r 4.3.0 \
-profile docker \
--mode fastq \
--input metadata_pcr_umi_airr_300.tsv \
--cprimers 's3://ngi-igenomes/test-data/airrflow/pcr_umi/cprimers.fasta' \
--vprimers 's3://ngi-igenomes/test-data/airrflow/pcr_umi/vprimers.fasta' \
--library_generation_method specific_pcr_umi \
--cprimer_position R1 \
--umi_length 15 \
--umi_start 0 \
--umi_position R1 \
-c resource.config \
--outdir bulk_fastq_results \
-resume
```

Of course you can wrap all your code in a bash file. We prepared one for you and it's available [here](https://github.com/nf-core/airrflow/tree/dev/docs/usage/bulk_tutorial/bulk_sample_code/airrflow_bulk_b_fastq.sh).
With the bash file, it's easy to run the pipeline with a single-line command.

```bash
bash airrflow_bulk_b_fastq.sh
```

If no UMI barcodes were used, set the `--library_generation_method specific_pcr`, and the UMI length will be set automatically to 0.

> [!TIP]
> Please ensure you modify the parameters when running the pipeline on your own data to match the specific details of your library preparation protocol.

> [!WARNING]
> When launching a Nextflow pipeline with the `-resume` option, any processes that have already been run with the exact same code, settings and inputs will be cached and the pipeline will resume from the last step that changed or failed with an error. The benefit of using "resume" is to avoid duplicating previous work and save time when re-running a pipeline.
> We include "resume" in our Nextflow command as a precaution in case anything goes wrong during execution. After fixing the issue, you can relaunch the pipeline with the same command, it will resume running from the point of failure, significantly reducing runtime and resource usage.

After launching the pipeline the following will be printed to the console output, followed by some the default parameters used by the pipeline and execution log of airrflow processes:

```bash
 N E X T F L O W   ~  version 24.10.5

WARN: It appears you have never run this project before -- Option `-resume` is ignored
Launching `https://github.com/nf-core/airrflow` [fabulous_cantor] DSL2 - revision: d91dd840f4 [4.3.0]


------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/airrflow 4.3.0
------------------------------------------------------
```

Once the pipeline has finished successfully, the following message will appear:

```bash
-[nf-core/airrflow] Pipeline completed successfully-
Completed at: 11-Mar-2025 15:01:11
Duration    : 23m 46s
CPU hours   : 5.9
Succeeded   : 271
```

## Understanding the results

After running the pipeline, several sub-folders are available under the results folder.

```bash
Airrflow_report.html
- fastp
- fastqc
- vdj_annotation
- qc_filtering
- clonal_analysis
- repertoire_comparison
- multiqc
- parsed_logs
- report_file_size
- pipeline_info
```

The summary report, named `Airrflow_report.html`, provides an overview of the analysis results, such as an overview of the number of sequences per sample in each of the pipeline steps, the V(D)J gene assignment and QC, and V gene family usage. Additionally, it contains links to detailed reports for other specific analysis steps.

The analysis steps and their corresponding folders, where the results are stored, are briefly listed below. Detailed documentation on the pipeline output can be found on the [Output documentation page](https://nf-co.re/airrflow/docs/output/).

1. Quality control
   - `fastp` is used to perform quality control, adapter trimming, quality filtering, per-read quality pruning of the FASTQ data. The results are stored under the folder 'fastp'.
   - `FastQC` is applied to do some quality control checks on raw sequence data. The fastqc report for each fastq file is under the folder named 'fastqc'. The aggregated QC reports for all samples can also be checked in the MultiQC report.

2. Sequence assembly
   - [`pRESTO`](https://presto.readthedocs.io/en/stable/) is a tool part of Immcantation for processing raw reads from high-throughput sequencing of B cell and T cell repertoires. It includes features for quality control, primer masking, annotation of reads with sequence embedded barcodes, generation of unique molecular identifier (UMI) consensus sequences, assembly of paired-end reads and identification of duplicate sequences. The 'presto' folder contains the intermediate pRESTO steps logs.

3. V(D)J annotation and filtering
   - In this step, V(D)J gene segments are inferred using the provided germline reference and [`IgBLAST`](https://www.ncbi.nlm.nih.gov/igblast/). Alignments are annotated in AIRR format. Non-productive sequences and sequences with low alignment quality are filtered out unless otherwise specified. The intermediate results are stored under the folder named 'vdj_annotation'.

4. Post alignment QC and filtering
   - Duplicates detected after alignment are collapsed in this step and the results are available under folder 'qc-filtering'.

5. Clonal analysis
   - Results of the clonal threshold determination using `SHazaM` should be inspected in the html report under the 'clonal_analysis/find_threshold' folder. If the automatic threshold is unsatisfactory, you can set the threshold manually and re-run the pipeline.
     (Tip: use -resume whenever running the Nextflow pipeline to avoid duplicating previous work).
   - Clonal inference is performed with `SCOPer`. Clonal inference results as well as clonal abundance and diversity plots can be inspected in the html report in the folder 'clonal_analysis/clonal_assignment'. For BCR sequencing data, mutation frequency is also computed using `SHazaM` at this step and plotted in the report. The `repertoires` subfolder contains the AIRR formatted files with the clonal assignments in a new column `clone_id` and mutation frequency in the column `mu_freq`. The `tables` subfolder contains the tabulated abundance and diversity computation as well as a table with the number of clones and their size. The `ggplots` subfolder contains the abundance and diversity plots as an `RData` object for loading and customization in R.
   - If lineage trees were computed using `Dowser`, a folder under 'clonal_analysis/dowser_lineages' will be present. The trees can be inspected in the html report and saved as PDF. Additionally, an `RDS` object with the formatted trees can also be loaded in R for customizing the lineage tree plots with Dowser.

6. Repertoire analysis
   - Comparison of several repertoire characteristics, such as V gene usage, across subjects, time points or cell populations. All associated plots and tables are available under the folder `repertoire_comparison`. The plots are also included in the `Airrflow_report.html` file. This report is generated from an R markdown `Rmd` file. It is possible to customize this to meet the user's needs by editing the report and then providing the edited Rmd file with the `--report_rmd` parameter. Check the remaining [Report parameters](https://nf-co.re/airrflow/parameters/#report-options) for further customizing the report.

7. Other reporting
   Additional reports are also generated, including:
   - `MultiQC` report: summarizes QC metrics across all samples.
   - `Pipeline_info` report: various reports relevant to the running and execution of the pipeline.
   - `Report_file_size` report: Summary of the number of sequences left after each of the most important pipeline steps.
   - `parsed_logs` report: Summary of the number of sequences left after each of the most important pipeline steps.

## Find out more

<<<<<<< HEAD
nf-core/airrflow is a standardized pipeline that performs the different computational analysis steps and provides standard figures for a first data exploration. The computations results (e.g. clonal inference, mutation frequency analysis) are stored in the output AIRR rearrangement repertoire files in newly generated columns under `clonal_analysis/clonal_assignment/all_repertoires`. You can use these Airrflow results as input for customized analyses using R and the Immcantation tools. You can find the tutorial for Immcantation's single-cell V(D)J analysis [here](https://immcantation.readthedocs.io/en/stable/getting_started/10x_tutorial.html).
=======
To continue learning about how to use nf-core/airrflow please check out the following documentation:
>>>>>>> 1daa0857bddc3046009a1002828b323edd565155

- [nf-core/airrflow usage documentation](https://nf-co.re/airrflow/docs/usage)
- [nf-core/airrflow parameters documentation](https://nf-co.re/airrflow/parameters)
- [FAQ page](./FAQ.md)

The nf-core troubleshooting documentation will also help you troubleshoot your Nextflow errors.

- [nf-core troubleshooting](https://nf-co.re/docs/usage/troubleshooting/overview)
