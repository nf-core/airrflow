# nf-core/bcellmagic: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [pRESTO](#presto) - read pre-processing
  * [Filter by sequence quality](#filter-by-sequence-quality) - filter sequences by quality
  * [Mask primers](#mask-primers) - Masking primers
  * [Pair mates](#pair-mates) - Pairing sequence mates.
  * [Cluster sets](#cluster-sets) - Cluster sequences according to similarity.
  * [Build consensus](#build-UMI-consensus) - Build UMI consensus
  * [Re-pair mates](#re-pair-mates) - Re-pairing sequence mates.
  * [Assemble mates](#assemble-mates) - Assemble sequence mates.
  * [Remove duplicates](#remove-duplicates) - Remove read duplicates.
  * [Filter sequences for at least 2 representative](#filter-sequences-for-at-least-2-representative) Filter sequences that do not have at least 2 reads assigned.
* [Change-O](#change-o) - Assign genes and clonotyping
  * [Assign genes with Igblast](#assign-genes-with-igblast)
  * [Make database from assigned genes](#make-database-from-assigned-genes)
  * [Removal of non-productive sequences](#removal-of-non-productive-sequences)
  * [Selection of IGH / TR sequences](#selection-of-IGH-/-TR-sequences)
  * [Convert database to fasta](#convert-database-to-fasta)
* [Shazam](#shazam) - Genotyping and Clonal threshold
  * [Genotyping and hamming distance threshold](#determining-genotype-and-hamming-distance-threshold)
* [Change-O define clones ](#change-o-define-clones)
  * [Define clones](#define-clones) - Defining clonal B-cell or T-cell groups
  * [Reconstruct germlines](#reconstruct-germlines) - Reconstruct gene calls of germline sequences
* [Lineage reconstruction](#lineage-reconstruction) - Clonal lineage reconstruction.
* [Repertoire analysis](#repertoire-analysis) - Repertoire analysis and comparison.
* [Log parsing](#log-parsing) - Log parsing.
* [Databases](#databases)
* [MultiQC](#MultiQC) - MultiQC
* [Pipeline information](#pipeline-information) - Pipeline information

## FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

**NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output files:**

* `fastqc/`
  * `*_fastqc.html`: FastQC report containing quality metrics for your untrimmed raw fastq files.
* `fastqc/zips/`
  * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

## presto

### Filter by sequence quality

Filters reads that are below a quality threshold by using the tool [FilterSeq](https://presto.readthedocs.io/en/version-0.5.11/tools/FilterSeq.html) from the pRESTO Immcantation toolset. The default quality threshold is 20.

**Output directory: `results/presto/01-filterseq`**

For each analyzed sample there is a subfolder containing:

* `logs`
  * Raw command logs of the process that will be parsed to generate a report.
* `tabs`
  * Table containing read ID and quality for each of the read files.

### Mask primers

Masks primers that are provided in the C-primers and V-primers input files. It uses the tool [MaskPrimers](https://presto.readthedocs.io/en/version-0.5.11/tools/MaskPrimers.html) of the pRESTO Immcantation toolset.

**Output directory: `results/presto/02-maskprimers`**

For each analyzed sample there is a subfolder containing:

* `logs`
  * Raw command logs of the process that will be parsed to generate a report.
* `tabs`
  * Table containing a read ID, the identified matched primer and the error for primer alignment.

### Pair mates

Pair read mates using [PairSeq](https://presto.readthedocs.io/en/version-0.5.11/tools/PairSeq.html) from the pRESTO Immcantation toolset.

**Output directory: `results/presto/03-pairseq`**

For each analyzed sample there is a subfolder containing:

* `logs`
  * Raw command logs of the process that will be parsed to generate a report.

### Cluster sets

Cluster sequences according to similarity, using [ClusterSets set](https://presto.readthedocs.io/en/version-0.5.11/tools/ClusterSets.html#clustersets-set). This step is introduced to deal with too low UMI diversity.

**Output directory: `results/presto/04-cluster_sets`**

For each analyzed sample there is a subfolder containing:

* `logs`
  * Raw command logs of the process that will be parsed to generate a report.
* `tabs`
  * Table containing a read ID, the identified barcode, the cluster id and the number of sequences in the cluster.

### Parse clusters

Annotate cluster ID as part of the barcode, using [Parseheaders copy](https://presto.readthedocs.io/en/stable/tools/ParseHeaders.html#parseheaders-copy). This step is introduced to deal with too low UMI diversity.

**Output directory: `results/presto/05-parse_clusters`**

For each analyzed sample there is a subfolder containing:

* `logs`
  * Raw command logs of the process that will be parsed to generate a report.

### Build UMI consensus

Build sequence consensus from all sequences that were annotated to have the same UMI. Uses [BuildConsensus](https://presto.readthedocs.io/en/version-0.5.11/tools/BuildConsensus.html) from the pRESTO Immcantation toolset.

**Output directory: `results/presto/06-build_consensus`**

For each analyzed sample there is a subfolder containing:

* `sample_command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `*.tab`
  * Table containing the sequence barcode, number of sequences used to build the consensus (SEQCOUNT), the identified primer (PRIMER), the number of sequences for each primer (PRCOUNT), the primer consensus (PRCONS), the primer frequency (PRFREQ) and the number of sequences used to build the consensus (CONSCOUNT).

### Re-pair mates

Re-pair read mates using [PairSeq](https://presto.readthedocs.io/en/version-0.5.11/tools/PairSeq.html) from the pRESTO Immcantation toolset.

**Output directory: `results/presto/07-pairseq_postconsensus`**

For each analyzed sample there is a subfolder containing:

* `logs`
  * Log of the process that will be parsed to generate a report.

### Assemble mates

Assemble read mates using [AssemblePairs](https://presto.readthedocs.io/en/version-0.5.11/tools/AssemblePairs.html) from the pRESTO Immcantation toolset.

**Output directory: `results/presto/08-assemblepairs`**

For each analyzed sample there is a subfolder containing:

* `logs`
  * Log of the process that will be parsed to generate a report.
* `tabs`
  * Parsed log contaning the sequence barcodes and sequence length, bases of the overlap, error of the overlap and p-value.

### Remove duplicates

Remove duplicates using [CollapseSeq](https://presto.readthedocs.io/en/version-0.5.11/tools/CollapseSeq.html) from the pRESTO Immcantation toolset.

**Output directory: `results/presto/09-collapseseq`**

For each analyzed sample there is a subfolder containing:

* `logs`
  * Log of the process that will be parsed to generate a report.
* `tabs`
  * Parsed log containing the sequence barcodes, header information and deduplicate count.

### Filter sequences for at least 2 representatives

Remove sequences which do not have 2 representative using [SplitSeq](https://presto.readthedocs.io/en/version-0.5.11/tools/SplitSeq.html) from the pRESTO Immcantation toolset.

**Output directory: `results/presto/10-splitseq`**

For each analyzed sample there is a subfolder containing:

* `logs`
  * Log of the process that will be parsed to generate a report.

## Change-O

### Assign genes with Igblast

Assign genes with Igblast, using the IMGT database is performed by the [AssignGenes](https://changeo.readthedocs.io/en/version-0.4.5/examples/igblast.html#running-igblast) command of the Change-O tool from the Immcantation Framework.

**Output directory: `results/changeo/01-assigngenes`**

For each analyzed sample there is a subfolder containing:

* `fasta/*.fasta`
  * Igblast results converted to fasta format with genotype V-call annotated in the header.

### Make database from assigned genes

A table is generated with [MakeDB](https://changeo.readthedocs.io/en/version-0.4.5/examples/igblast.html#processing-the-output-of-igblast).

**Output directory: `results/changeo/02-makedb`**

For each analyzed sample there is a subfolder containing:

* `logs`
  * Log of the process that will be parsed to generate a report.
* `tab`
  * Table in AIRR format containing the assigned gene information and metadata provided in the starting metadata sheet.

### Removal of non-productive sequences

Non-functional sequences are removed with [ParseDb](https://changeo.readthedocs.io/en/version-0.4.5/tools/ParseDb.html).

**Output directory: `results/changeo/03-parsedb_split`**

For each analyzed sample there is a subfolder containing:

* `logs`
  * Log of the process that will be parsed to generate a report.
* `tab`
  * Table in AIRR format containing the assigned gene information, with only productive sequences and metadata provided in the starting metadata sheet.

### Selection of IGH / TR sequences

Heavy chain sequences (IGH) are selected if 'ig' locus is selected, TR sequences are selected if 'tr' locus is selected. The tool [ParseDb](https://changeo.readthedocs.io/en/version-0.4.5/tools/ParseDb.html) is employed.

**Output directory: `results/changeo/04-parsedb_select`**

For each analyzed sample there is a subfolder containing:

* `logs`
  * Log of the process that will be parsed to generate a report.
* `tab`
  * Table in AIRR format containing the assigned gene information, with only productive sequences and IGH/TR sequences, and metadata provided in the starting metadata sheet.

### Convert database to fasta

Sequences in are additionally converted to a fasta file with the [ConvertDb](https://changeo.readthedocs.io/en/version-0.4.5/tools/ConvertDb.html?highlight=convertdb) tool.

**Output directory: `results/changeo/05-convertdb-fasta`**

For each analyzed sample there is a subfolder containing:

* `fasta`
  * Fasta file containing the processed sequences with the barcode ID and allele annotation in the header.

## Shazam

### Determining genotype and hamming distance threshold

Determining genotype and the hamming distance threshold of the junction regions for clonal determination using the [tigGER](https://tigger.readthedocs.io/en/0.3.1/) and [Shazam](https://shazam.readthedocs.io/en/version-0.1.11_a/).

**Output directory: `results/shazam/02-genotyping`**

For each subject (patient) there is a subfolder containing:

* `threshold`
  * Hamming distance threshold of the Junction regions as determined by Shazam.
* `plots`
  * Plot of the Hamming distance distribution between junction regions displaying the threshold for clonal assignment as determined by Shazam.
* `genotype`
  * `genotype.pdf`: Plot representing the patient genotype assessed by TigGER.
  * `v_genotype.fasta`: Fasta sequences of the personalized patient genotype.

## Change-O define clones

### Define clones

Assigning clones to the sequences obtained from IgBlast with the [DefineClones](https://changeo.readthedocs.io/en/version-0.4.5/tools/DefineClones.html?highlight=DefineClones) Immcantation tool.

**Output directory: `results/changeo/06-define_clones`**

For each subject (patient) there is a subfolder containing:

* `tab`
  * Table in AIRR format containing the assigned gene information and an additional field with the clone id.

### Reconstruct germlines

Reconstructing the germline sequences with the [CreateGermlines](https://changeo.readthedocs.io/en/version-0.4.5/tools/CreateGermlines.html#creategermlines) Immcantation tool.

**Output directory: `results/changeo/07-create-germlines`**

For each subject (patient) there is a subfolder containing:

* `tab`
  * Table in AIRR format contaning the assigned gene information and an additional field with the germline reconstructed gene calls.

## Lineage reconstruction

Reconstructing clonal linage with the [Alakazam R package](https://alakazam.readthedocs.io/en/stable/) from the Immcantation toolset.

**Output directory: `results/lineage_reconstruction`**

For each subject (patient) there is a subfolder containing the following:

* `tab`
  * `Clones_table_patient.tsv`: contains a summary of the clones found for the patient, and the number of unique and total sequences identified in each clone.
  * `Clones_table_patient_filtered_between_3_and_1000.tsv`: contains a summary of the clones found for the patient, and the number of unique and total sequences identified in each clone, filtered by clones of size between 3 and 1000, for which the lineages were reconstructed and the trees plotted.
  * `xxx_germ-pass.tsv`: AIRR format table with all the sequences from a patient after the germline annotation step.
* `Clone_tree_plots`
  * Contain a rooted graphical representation of each of the clones, saved in pdf format.
* `Graphml_trees`
  * All lineage trees for the patient exported in a GraphML format: `All_graphs_patient.graphml`.

## Repertoire comparison

Calculation of several repertoire characteristics (diversity, abundance, V gene usage) for comparison between subjects, time points and cell populations. An Rmarkdown report is generated with the [Alakazam R package](https://alakazam.readthedocs.io/en/stable/).

**Output directory: `results/repertoire_analysis`**

* `repertoire_comparison/`
  * `all_data.tsv`: AIRR format table containing the processed sequence information for all subjects.
  * `Abundance`: contains clonal abundance calculation plots and tables.
  * `Diversity`: contains diversity calculation plots and tables.
  * `V_family`: contains V gene and family distribution calculation plots and tables.
* `Bcellmagic_report.html`
  * Contains the repertoire comparison results in an html report form: Abundance, Diversity, V gene usage tables and plots. Comparison between treatments and subjects.

## Log parsing

Parsing the logs from the previous processes.

**Output directory: `results/parsed_logs`**

* `sequences_table`: table summarizing of the number of sequences after the most important pipeline steps.

## Databases

Copy of the downloaded IMGT database by the process `fetch_databases`, used for the gene assignment step.

## MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability.

For more information about how to use MultiQC reports, see [https://multiqc.info](https://multiqc.info).

**Output files:**

* `multiqc/`
  * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  * `multiqc_plots/`: directory containing static images from the report in various formats.

## Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files:**

* `pipeline_info/`
  * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  * Documentation for interpretation of results in HTML format: `results_description.html`.
