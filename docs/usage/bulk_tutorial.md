# nf-core/airrflow: Bulk AIRRseq tutorial

This tutorial provides a step by step introduction on how to run nf-core/airrflow on bulk AIRR-seq data.

## Pre-requisites

> [!INSTALLATION]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set up Nextflow and a container engine needed to run this pipeline. At the moment, nf-core/airrflow does NOT support using conda virtual environments for dependency management, only containers are supported. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) before running the workflow on actual data.

For the purpose of running this tutorial on your local machine, we recommend a docker installation.

To install docker, follow the instructions [here](https://docs.docker.com/engine/install/). After installation docker on linux, don't forget to check the [post-installation steps](https://docs.docker.com/engine/install/linux-postinstall/).

## Testing the pipeline with built-in tests

Once you have set up your Nextflow and container (docker or singularity), test the airrflow pipeline with built-in test.

```bash
nextflow run nf-core/airrflow -r 4.2.0 -profile test,docker --outdir test_results
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
- A configuration file specifying the system's maximum available RAM memory, CPUs and running time. This will ensure that no pipeline process requests more resources than available in the compute infrastructure where the pipeline is running. The resource configuration file is provided with the `-c` option. In this example we set the maximum RAM memory to 20GB, we restrict the pipeline to use 8 CPUs and to run for a maximum of 24 hours. Depending on the size of your dataset, it might be required to extend the running time. You can also remove the "time" parameter from the configuration file to allow for unlimited runtime.

```bash
process {
   resourceLimits = [cpus: 8, memory: 20.GB, time: 24.h]
}
```

> [Tip]
> Before setting memory and cpus in the configuration file, we recommend verifying the available memory and cpus on your system. Otherwise, exceeding the system's capacity may result in an error indicating that you requested more cpus than available or run out of memory.

> [Tip]
> When running nf-core/airrflow with your own data, provide the full path to your input files under the filename column.

A prepared samplesheet for this tutorial can be found [here](bulk_tutorial/bulk_sample_code/metadata_pcr_umi_airr_300.tsv), and the configuration file is available [here](bulk_tutorial/bulk_sample_code/resource.config).
Download both files to the directory where you intend to run the airrflow pipeline.

## Choosing the right protocol profile

Bulk BCR and TCR targeted sequencing can be performed with a wide variety of protocols, using different library preparation methods. Different protocols usually use different amplification primers, UMI barcode lengths and position, which require different parameter setting to run Airrflow pipeline. To make it easier to run the pipeline on commonly used commercially available kits, we provide parameter presets as profiles. A full [list of protocol profiles](https://nf-co.re/airrflow/docs/usage/#supported-protocol-profiles) is available on the usage documentation page.

You can provide a protocol profile with the `-profile` parameter, followed by other profiles, such as the container engine profile in a comma separated fashion. You will then usually only need to provide as additional parameters the input samplesheet, resource config file and output directory path. However, if you want to override any option or add additional parameters, you can provide them to the airrflow launching command as any parameters in the launch command will override the parameters in the profile.

```bash
nextflow run nf-core/airrflow -r 4.2.0 \
-profile <protocol-profile-name>,docker \
--input samplesheet.tsv \
-c resource.config \
--outdir bulk_fastq_results \
-resume
```

> [Tip]
> We're always looking forward to expanding the set of protocol profiles readily available for other users. Feel free to open an issue and create a pull request to add a new profile that you want to share with other users or ask in the nf-core `#airrflow` [slack channel](https://nf-co.re/join) if you have any questions in doing so.

## Analyzing a dataset with a custom library preparation method

If your dataset was generated using a custom library preparation method, you can manually set the relevant parameters according to your protocol design, similar to the approach we used for the samples in this tutorial. For more examples on how to set the parameters for custom protocols check the [usage documentation](https://nf-co.re/airrflow/docs/usage/#supported-bulk-library-generation-methods-protocols) page.

The BCRseq dataset used in this tutorial was obtained with a multiplexed PCR protocol using custom C-region and V-region primers. We stored the sequences for the V-region primers as well as the C-region primers in AWS S3, and the links are provided in the Nextflow command which will be fetched by nextflow automatically when executing the command. You can also provide the full path to the custom primers fasta files.

The command to launch nf-core/airrflow for the dataset in this tutorial is the following:

```bash
nextflow run nf-core/airrflow -r 4.2.0 \
-profile docker \
--mode fastq \
--input samplesheet.tsv \
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

Of course you can wrap all your code in a bash file. We prepared one for you and it's available [here](bulk_tutorial/bulk_sample_code/airrflow_bulk_b_fastq.sh).
With the bash file, it's easy to run the pipeline with a single-line command.

```bash
bash airrflow_bulk_b_fastq.sh
```

If no UMI barcodes were used, set the `--library_generation_method specific_pcr`, and the UMI length will be set automatically to 0.

> [Warning!]
> Please ensure you modify the parameters when running the pipeline on your own data to match the specific details of your library preparation protocol.

> [Tip]
> When launching a Nextflow pipeline with the `-resume` option, any processes that have already been run with the exact same code, settings and inputs will be cached and the pipeline will resume from the last step that changed or failed with an error. The benefit of using "resume" is to avoid duplicating previous work and save time when re-running a pipeline.
> We include "resume" in our Nextflow command as a precaution in case anything goes wrong during execution. After fixing the issue, you can relaunch the pipeline with the same command, it will resume running from the point of failure, significantly reducing runtime and resource usage.

After launching the pipeline the following will be printed to the console output, followed by some Nextflow parameters and executions of Airrflow processes:

```bash
 N E X T F L O W   ~  version 24.10.5

WARN: It appears you have never run this project before -- Option `-resume` is ignored
Launching `https://github.com/nf-core/airrflow` [fabulous_cantor] DSL2 - revision: d91dd840f4 [4.2.0]


------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/airrflow 4.2.0
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

## Important considerations for clonal analysis

An important step in the analysis of AIRR sequencing data is inferring B cell and T cell clones, or clonal groups, sometimes also called clonotypes. These are cells that are derived from the same progenitor cell through clonal expansion. For T cells, this definition is more strict as T cells do not undergo somatic hypermutation, so the TCRs from T cells in the same clone should be identical. For B cells, on the other hand, the BCRs from cells in the same clone can differ due to somatic hypermutation. They also can have a variety of isotypes.

There are two crucial considerations when defining clonal groups with nf-core/airrflow: across which samples should clonal groups be defined, and what should be the clonal threshold, i.e. how different can these receptors be, so that these are assigned to the same clonal group. These are discussed in detail in the following sections.

### Defining clonal groups across samples

Often times we want to analyze clonal groups from the same individual or animal model across time, different conditions or across samples extracted from different tissues. To ensure that the same clone ID (field `clone_id` in the output AIRR rearrangement file) is assigned to the same BCR / TCR clone across these conditions to be able to track the clones, the clonal inference step should be done pulling the sequences from these samples together. This is why, by default, nf-core/airrflow uses the `subject_id` column to group samples prior to defining clonal groups, so it is important to set the exact same subject ID to samples from the same individual across different conditions.

The sample grouping can also be controlled with the [`--cloneby`](https://nf-co.re/airrflow/4.2.0/parameters/#cloneby) parameter, by providing the name of the column containing the group information that should be used to pull the samples together before defining clonal groups (samples or rows with the same string in this column will be grouped together). You can create a new column if you wish for this purpose.

### Clonal inference method

nf-core/airrflow utilizes the Hierarchical clustering method in the [SCOPer](https://scoper.readthedocs.io/) Immcantation tool to infer clonal groups, which initially partitions the BCR / TCR sequences according to V gene, J gene and junction length. Then, it defines clonal groups within each partition by performing hierarchical clustering of the sequences within a partition and cutting the clusters according to an automatically detected or user-defined threshold. More details about this method can be found on the respective SCOPer [vignette](https://scoper.readthedocs.io/en/stable/vignettes/Scoper-Vignette/#). Details on how to determine the clonal threshold can be found in the next section.

### Setting a clonal threshold

The clonal threshold can also be customized through the `--clonal_threshold` parameter. By default, `--clonal_threshold` is set to be 'auto', allowing the threshold of how different two BCRs - or specifically their junction regions - can be to be assigned to the same clonal to be determined automatically using a method included in the [SHazaM](https://shazam.readthedocs.io/) Immcantation tool. You can read more details about the method in the SHazaM [vignette](https://shazam.readthedocs.io/en/stable/vignettes/DistToNearest-Vignette/).

For BCR data, we recommend using this default setting initially. After running the pipeline, you can review the automatically calculated threshold in the `find_threshold` report to make sure it is fitting the data appropriately. If the threshold is unsatisfactory, you can re-run the pipeline with a manually specified threshold (e.g. `--clonal_threshold 0.1`) that is appropriate for your data. For a low number of sequences that are insufficient to satisfactorily determine a threshold with this method, we generally recommend a threshold of 0.1 (length-normalized Hamming distance of nearest neighbors) for human BCR data.

Since TCRs do not undergo somatic hypermutation, TCR clones are defined strictly by identical junction regions. For this reason, the `--clonal_threshold` parameter should be set to 0 for TCR data.

### Including BCR lineage tree computation

BCR lineage tree computation is performed using the [Dowser](https://dowser.readthedocs.io/) Immcantation package. This step is skipped by default because it can be time-consuming depending on the size of the input data and the size of the clonal groups. To enable lineage tree computation, add the `--lineage_trees` parameter set to true. You can easily add lineage tree computation to a previous analysis by re-running the pipeline with the `-resume` so all the previous analysis steps are cached and not recomputed.

Dowser supports different methods for the lineage tree computation, `raxml` is the default but you can set other methods with the `--lineage_tree_builder` parameter, and provide the software executable with the `--lineage_tree_exec` parameter.

## Understanding the results

After running the pipeline, several subfolders are available under the results folder.

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
   - Clonal inference is performed with `SCOPer`. Clonal inference results as well as clonal abundance and diversity plots can be inspected in the html report in the folder 'clonal_analysis/define_clones'. For BCR sequencing data, mutation frequency is also computed using `SHazaM` at this step and plotted in the report. The `repertoires` subfolder contains the AIRR formatted files with the clonal assignments in a new column `clone_id` and mutation frequency in the column `mu_freq`. The `tables` subfolder contains the tabulated abundance and diversity computation as well as a table with the number of clones and their size. The `ggplots` subfolder contains the abundance and diversity plots as an `RData` object for loading and customization in R.
   - If lineage trees were computed using `Dowser`, a folder under 'clonal_analysis/dowser_lineages' will be present. The trees can be inspected in the html report and saved as PDF. Additionally, an `RDS` object with the formatted trees can also be loaded in R for customizing the lineage tree plots with Dowser.

6. Repertoire analysis

   - Example calculation of several repertoire characteristics, e.g. V gene usage, for comparison between subjects, time points or cell populations is shown in the html report under `repertoire_comparison`. This report is generated from an Rmarkdown `Rmd` file. It is possible to customize this to fit the user's needs by editing the report and then providing the edited Rmd file with the `--report_rmd` parameter. Check also the remaining [Report parameters](https://nf-co.re/airrflow/parameters/#report-options) for further customizing the report.

7. Other reporting
   Additional reports are also generated, including:
   - `MultiQC` report: summarizes QC metrics across all samples.
   - `Pipeline_info` report: various reports relevant to the running and execution of the pipeline.
   - `Report_file_size` report: Summary of the number of sequences left after each of the most important pipeline steps.
   - `parsed_logs` report: Summary of the number of sequences left after each of the most important pipeline steps.

## Costumizing your analysis and generating your own figures

nf-core/airrflow is a standardized pipeline that performs the different computational analysis steps and provides standard figures for a first data exploration. The computations results (e.g. clonal inference, mutation frequency analysis) are stored in the output AIRR rearrangement repertoire files in newly generated columns under `clonal_analysis/define_clones/all_repertoires`. You can use these Airrflow results as input for customized analyses using R and the Immcantation tools. You can find the tutorial for Immcantation's single-cell V(D)J analysis [here](https://immcantation.readthedocs.io/en/stable/getting_started/10x_tutorial.html).

## Updating process resource requests

By default the pipeline has set reasonable process resource requests (number of CPUs, RAM memory, time limits) to the compute system. Depending on the size of your datasets or your running infrastructure you can customize these requests to meet your needs.

To update the resource requests for a specific pipeline process, you can do so in the `resource.config` file provided with the `-c` parameter. For example, to update the resource requests for the `CHANGEO_ASSIGNGENES` process:

````bash
process {
   resourceLimits = [cpus: 8, memory: 72.GB, time: 24.h]

   withName:CHANGEO_ASSIGNGENES {
        cpus   = 2
        memory = 10.GB
        time   = 5h
   }
}

In nf-core pipelines, each process has a label indicating the resources that are being requested (`process_low`, `process_medium`, `process_high`, ...). The CPUs, RAM and time set up for each of these labels can be found in the [base.config](../../../conf/base.config) file. You can update the resource requests for all processes with a specific label by adding a new setting in your `resource.config` file provided with the `-c` parameter. For example here we update the process requests of processes with the `process_high` label:

```bash
process {
   resourceLimits = [cpus: 24, memory: 100.GB, time: 24.h]

   withLabel:process_high {
        cpus   = 24
        memory = 100.GB
        time   = 10h
   }
}
````

Note that the resource requests will never exceed what is specified in the `resourceLimits` line, so if you do want to increase the resource requests for specific processes, you should also increase the `resourceLimits` requests and run the pipeline in a compute infrastructure with sufficient resources. In this exmaple we also have updated the `resourceLimits` to reflect that.

> [Tip]
> For more information about nf-core pipeline resource configurations, check out the [nf-core pipeline configuration docs](https://nf-co.re/docs/usage/getting_started/configuration).
