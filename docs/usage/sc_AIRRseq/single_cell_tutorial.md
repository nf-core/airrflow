# nf-core/airrflow: Single-cell tutorial

This tutorial provides a step by step introduction on how to run nf-core/airrflow on single-cell BCR-seq data or single-cell TCR-seq data.

## Pre-requisites

> [!INSTALLATION]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set up Nextflow and a container engine needed to run this pipeline. At the moment, nf-core/airrflow does NOT support using conda virtual environments for dependency management, only containers are supported. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) before running the workflow on actual data.

For the purpose of running this tutorial on your local machine, we recommend a docker installation.

To install docker, follow the instructions [here](https://docs.docker.com/engine/install/). After docker installation on linux system, don't forget to check the [post-installation steps](https://docs.docker.com/engine/install/linux-postinstall/).

## Testing the pipeline with built-in tests

Once you have set up your Nextflow and container (docker or singularity), test the airrflow pipeline with built-in test. 

```bash
nextflow run nf-core/airrflow -r 4.2.0 -profile test,docker --outdir test_results
```


## Running airrflow pipeline from two different input formats
There are two acceptable input formats for airrflow single-cell AIRRseq pipeline: AIRR rearrangement or fastq format. 

For this tutorial we will practice on both of the input formats. 

## AIRR rearrangement format
### Datasets

For this tutorial we will use subsampled PBMC single-cell BCR sequencing data from two subjects, before (d0) and after flu vaccination (d12). The dataset is available on [Zenodo](https://zenodo.org/doi/10.5281/zenodo.11373740). You don't need to download the samples bacause the links to the samples are already in the samplesheet we provide. 

### Preparing samplesheet and configuration file

To run the pipeline, a tab-seperated samplesheet and a configuration file must be prepared. Details on the requied columns of a samplesheet are available [here](https://nf-co.re/airrflow/usage#input-samplesheet). 

A prepared samplesheet for this tutorial can be found [here](sample_data_code/assembled_samplesheet.tsv), and the configuration file is available [here](sample_data_code/resource.config). 
Download both files to the directory where you intend to run the airrflow pipeline. 

> [Tip] 
> Before setting memory and cpus in the configuration file, we recommend verifying the available memory and cpus on your system. Otherwise, exceeding the system's capacity may result in unexpected errors. 

### Running airrflow 

With all the files ready, you can proceed to run the airrflow pipeline. 

```bash
nextflow run nf-core/airrflow -r 4.2.0 \
-profile docker \
--mode assembled \
--input assembled_samplesheet.tsv \
--outdir sc_from_assembled_results  \
-c resource.config \
-resume
```
Of course you can wrap all your code in a bash file. We prepared one for you and it's available [here](sample_data_code/airrflow_sc_from_assembled.sh).
With the bash file, it's easy to run the pipeline with a single-line command. 

```bash
bash airrflow_sc_from_assembled.sh
```


> [Tip]
> When launching a Nextflow pipelien with -resume option, any processes that have already been run with the exact same code, settings and inputs will be skipped. The benefit of using -resume is to avoid duplicating previous work and save time when re-running a pipeline.
> We include -resume in our Nextflow command  as a precaution in case anything goes wrong during execution. After fixing the issue, you can relaunch the pipeline with the same command, it will resume running from the point of failure, significantly reducing runtime and resource usage.   


## Fastq format
### Datasets
For this tutorial we will use subsampled blood single-cell TCR sequencing data of one subject generated from the 10x Genomic platform. The links to the fastq files are in the samplesheet. 

### Preparing samplesheet, gene reference and configuration file
To run the airrflow pipeline on single cell TCR or BCR sequencing data from fastq files, we need to prepare samplesheet, gene reference and configuration file in advance. Details on the requied columns of a samplesheet are available [here](https://nf-co.re/airrflow/usage#input-samplesheet). 

The prepared samplesheet for this tutorial is [here](sample_data_code/10x_sc_raw.tsv) and a prepared configuration file is [here](sample_data_code/resource.config). Download these two files to the directory where you intend to run the airrflow pipeline.

> [Tip] 
> Before setting memory and cpus in the configuration file, we recommend verifying the available memory and cpus on your system. Otherwise, exceeding the system's capacity may result in unexpected errors. 

Gene reference can be accessed at the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/downloads). Both human and mouse V(D)J references are available. Download the reference that corresponds to the species of your dataset. 



### Running airrflow
With all the files ready, it's time to run the airrflow pipeline. 

```bash
nextflow run nf-core/airrflow -r 4.2.0 \
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

Of course you can wrap all your code in a bash file. We prepared one for you and it's available [here](sample_data_code/airrflow_sc_from_fastq.sh).
With the bash file, it's easy to run the pipeline with a single-line command. 

```bash
bash airrflow_sc_from_fastq.sh
```

By default, clonal_threshold is set to be 'auto', allowing the Hamming distance threshold of junction regions to be determined automatically. For BCR data, we recommend using this default setting. After running the pipeline, review the automatically calculated threshold to make sure it is appropriate. If the threshold is unsatisfactory, you can re-run the pipeline with a manually specified clonal_threshold. 
In this tutorial, since the samples are TCRs, which do not have somatic hypermutation, clones are defined strictly by identical CDR3s. For this reason, we set the clone-threshold parameter to 0. 


## Understanding the results

After running the pipeline, several reports are generated under the result folder. 

![example of result folder](tutorial_images/airrflow_result_folder_example.png)


The summary report, named 'Airrflow_report.html', provides an overview of the analysis results, such as the V(D)J gene assignment and QC, and V gene family usage. Additionally, it contains links to detailed reports for other specific analysis steps. 

The analysis steps and their corresponding folders, where the results are stored, are listed below. 

1. QC and sequence assembly (if starting from fastq files). 
   - In this first step, Cell Ranger's VDJ algorithm is employed to assemble contigs, annotate contigs, call cells and generate clonoytpes. The results are stored in the 'cellranger' folder.  

2. V(D)J annotation and filtering. 
   - In this step, gene segments are assigned using a germline reference. Alignments are annotated in AIRR format. Non-productive sequences and sequences with low alignment quality are removed. Metadata is added. The results are under the folder named 'vdj_annotation'. 

3. QC filtering. 
   - In this step, cells without heavy chains or with multiple heavy chains are removed. Sequences in different samples that share the same cell_id and necleotide sequence are filtered out. The result are stored in the 'qc-filtering' folder. 

4. Clonal analysis. 
   - In this step, the Hamming distance threshold of the junction regions is determined when clonal_threshold is set to 'auto' (by default).it should be reviewed for accuracy once the result is out. The threshold result can be found under the folder clonal_analysis/find_threshold. 
   - If the automatic threshold is unsatisfactory, you can set the threshold manually and re-run the pipeline. 
   (Tip: use -resume whenever running the Nextflow pipeline to avoid duplicating previous work). 
   - For TCR data, where somatic hypermutation does not occur, set the clonal_threshold to 0 when running the Airrflow pipeline.  
   - Once the threshold is established, clones are assigned to the sequences. A variety of tables and plots associated with clonal analysis were added to the folder 'clonal_analysis/define_clones', such as  sequences_per_locus_table, sequences_per_c_call_table, sequences_per_constant_region_table,num_clones_table, clone_sizes_table,clone size distribution plot, clonal abundance plot, diversity plot and etc. 

5. Repertoire analysis. 
   - The output folder is'repertoire_comparison'. V gene distribution tables and plots are included in this folder.

6. Other reporting.
   - Additional reports are also generated, including: a multiqc report which summarizes QC metrics across all samples, pipeline_info reports and report_file_size reports.



## Including lineage tree computation

Lineage tree computation is skipped by default because it's time-consuming. To enable lineage tree computation, re-run the pipeline with the --lineage_trees parameter set to true. Remember to include the -resume parameter to avoid duplicating previous work.


## Downstream analysis

Airrflow is a standardized pipeline that is not highly flexible for customized downstream analysis. For such cases, you can use the Airrflow results as input for customized analyses using the Immcantation packages with appropriate parameters. You can find the tutorial for Immcantation's single-cell V(D)J analysis [here](https://immcantation.readthedocs.io/en/stable/getting_started/10x_tutorial.html). 


