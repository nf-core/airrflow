# nf-core/airrflow: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

TODO: update this to add/remove lines

- [FastP](#fastp) - read quality control, adapter trimming and read clipping
- [pRESTO](#presto) - read pre-processing
  - [Filter by sequence quality](#filter-by-sequence-quality) - filter sequences by quality
  - [Mask primers](#mask-primers) - Masking primers
  - [Pair mates](#pair-mates) - Pairing sequence mates.
  - [Cluster sets](#cluster-sets) - Cluster sequences according to similarity.
  - [Build consensus](#build-UMI-consensus) - Build consensus of sequences with the same UMI barcode.
  - [Re-pair mates](#re-pair-mates) - Re-pairing sequence mates.
  - [Assemble mates](#assemble-mates) - Assemble sequence mates.
  - [Remove duplicates](#remove-duplicates) - Remove and annotate read duplicates.
  - [Filter sequences for at least 2 representative](#filter-sequences-for-at-least-2-representative) Filter sequences that do not have at least 2 duplicates.
- [FastQC](#fastqc) - read quality control post-assembly
- [Change-O](#change-o) - Assign genes and clonotyping
  - [Assign genes with Igblast](#assign-genes-with-igblast)
  - [Make database from assigned genes](#make-database-from-assigned-genes)
  - [Removal of non-productive sequences](#removal-of-non-productive-sequences)
  - [Selection of IGH / TR sequences](#selection-of-IGH-/-TR-sequences)
  - [Convert database to fasta](#convert-database-to-fasta)
- [Shazam](#shazam) - Genotyping and Clonal threshold
  - [Genotyping and hamming distance threshold](#determining-hamming-distance-threshold)
- [Change-O define clones](#change-o-define-clones)
  - [Define clones](#define-clones) - Defining clonal B-cell or T-cell groups
  - [Reconstruct germlines](#reconstruct-germlines) - Reconstruct gene calls of germline sequences
- [Lineage reconstruction](#lineage-reconstruction) - Clonal lineage reconstruction.
- [Repertoire analysis](#repertoire-analysis) - Repertoire analysis and comparison.
- [Log parsing](#log-parsing) - Log parsing.
- [Databases](#databases)
- [MultiQC](#MultiQC) - MultiQC
- [Pipeline information](#pipeline-information) - Pipeline information

## Fastp

<details markdown="1">
<summary>Output files</summary>

- `fastp/`
  - `<sample_id>/`
    - `*.fastp.html`: Fast report containing quality metrics for the mated and quality filtered reads.
    - `*.fastp.json`: Zip archive containing the FastQC report, tab-delimited data file and plot images for the mated and quality filtered reads.
    - `*.fastp.log`: Fastp

</details>

[fastp](https://doi.org/10.1093/bioinformatics/bty560) gives general quality metrics about your sequenced reads, as well as allows filtering reads by quality, trimming adapters and clipping reads at 5' or 3' ends. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [fastp documentation](https://github.com/OpenGene/fastp).

## presto

> **NB:** If using the sans-UMI subworkflow by specifying `umi_length=0`, the presto directory ordering numbers will differ e.g., mate pair assembly results will be output to `presto/01-assemblepairs/<sampleID>` as this will be the first presto step.

### Filter by sequence quality

<details markdown="1">
<summary>Output files</summary>

- `presto/01-filterseq/<sampleID>`
  - `logs`: Raw command logs of the process that will be parsed to generate a report.
  - `tabs`: Table containing read ID and quality for each of the read files.

</details>

Filters reads that are below a quality threshold by using the tool [FilterSeq](https://presto.readthedocs.io/en/stable/tools/FilterSeq.html) from the pRESTO Immcantation toolset. The default quality threshold is 20.

### Mask primers

<details markdown="1">
<summary>Output files</summary>

- `presto/02-maskprimers/<sampleID>`
  - `logs`: Raw command logs of the process that will be parsed to generate a report.
  - `tabs`: Table containing a read ID, the identified matched primer and the error for primer alignment.

</details>

Masks primers that are provided in the C-primers and V-primers input files. It uses the tool [MaskPrimers](https://presto.readthedocs.io/en/stable/tools/MaskPrimers.html) of the pRESTO Immcantation toolset.

### Pair mates

<details markdown="1">
<summary>Output files</summary>

- `presto/03-pairseq/<sampleID>`
  - `logs`: Raw command logs of the process that will be parsed to generate a report.

</details>

Pair read mates using [PairSeq](https://presto.readthedocs.io/en/stable/tools/PairSeq.html) from the pRESTO Immcantation toolset.

### Cluster sets

<details markdown="1">
<summary>Output files</summary>

- `presto/04-cluster_sets/<sampleID>`
  - `logs`: Raw command logs of the process that will be parsed to generate a report.
  - `tabs`: Table containing a read ID, the identified barcode, the cluster id and the number of sequences in the cluster.

</details>

Cluster sequences according to similarity, using [ClusterSets set](https://presto.readthedocs.io/en/stable/tools/ClusterSets.html#clustersets-set). This step is introduced to deal with too low UMI diversity.

### Parse clusters

<details markdown="1">
<summary>Output files</summary>

- `presto/05-parse_clusters/<sampleID>`
  - `logs`: Raw command logs of the process that will be parsed to generate a report.

</details>

Annotate cluster ID as part of the barcode, using [Parseheaders copy](https://presto.readthedocs.io/en/stable/tools/ParseHeaders.html#parseheaders-copy). This step is introduced to deal with too low UMI diversity.

### Build UMI consensus

<details markdown="1">
<summary>Output files</summary>

- `presto/06-build_consensus/<sampleID>`
  - `logs`: Raw command logs of the process that will be parsed to generate a report.
  - `tabs`: Table containing the sequence barcode, number of sequences used to build the consensus (SEQCOUNT), the identified primer (PRIMER), the number of sequences for each primer (PRCOUNT), the primer consensus (PRCONS), the primer frequency (PRFREQ) and the number of sequences used to build the consensus (CONSCOUNT).

</details>

Build sequence consensus from all sequences that were annotated to have the same UMI. Uses [BuildConsensus](https://presto.readthedocs.io/en/stable/tools/BuildConsensus.html) from the pRESTO Immcantation toolset.

### Re-pair mates

<details markdown="1">
<summary>Output files</summary>

- `presto/07-pairseq_postconsensus/<sampleID>`
  - `logs`: Raw command logs of the process that will be parsed to generate a report.

</details>

Re-pair read mates using [PairSeq](https://presto.readthedocs.io/en/stable/tools/PairSeq.html) from the pRESTO Immcantation toolset.

### Assemble mates

<details markdown="1">
<summary>Output files</summary>

- `presto/08-assemblepairs/<sampleID>`
  - `logs`: Raw command logs of the process that will be parsed to generate a report.
  - `tabs`: Parsed log contaning the sequence barcodes and sequence length, bases of the overlap, error of the overlap and p-value.

</details>

Assemble read mates using [AssemblePairs](https://presto.readthedocs.io/en/stable/tools/AssemblePairs.html) from the pRESTO Immcantation toolset.

### Remove duplicates

<details markdown="1">
<summary>Output files</summary>

- `presto/09-collapseseq/<sampleID>`
  - `logs`: Raw command logs of the process that will be parsed to generate a report.
  - `tabs`: Parsed log containing the sequence barcodes, header information and deduplicate count.

</details>

Remove duplicates using [CollapseSeq](https://presto.readthedocs.io/en/stable/tools/CollapseSeq.html) from the pRESTO Immcantation toolset.

### Filter sequences for at least 2 representatives

<details markdown="1">
<summary>Output files</summary>

- `presto/10-splitseq/<sampleID>`
  - `logs`: Raw command logs of the process that will be parsed to generate a report.

</details>

Remove sequences which do not have 2 representative using [SplitSeq](https://presto.readthedocs.io/en/stable/tools/SplitSeq.html) from the pRESTO Immcantation toolset.

## FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `postassembly/`
    - `*_ASSEMBLED_fastqc.html`: FastQC report containing quality metrics for the mated and quality filtered reads.
    - `*_ASSEMBLED_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images for the mated and quality filtered reads.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** Two sets of FastQC plots are displayed in the MultiQC report: first for the raw _untrimmed_ and unmated reads and secondly for the assembled and QC filtered reads (but before collapsing duplicates). They may contain adapter sequence and potentially regions with low quality.

## Change-O

### Convert input to fasta, if needed

<details markdown="1">
<summary>Output files. Optional. </summary>

- `vdj_annotation/convert-db/<sampleID>`
  - `*.fasta`: The sequences in fasta format.
  - `*log.txt`: Log of the process that will be parsed to generate a report.

</details>

This folder is generated when the input data are AIRR-C formatted rearrangement tables that need to
be reprocessed (`--reassign true`). For example, 10x Genomics' `airr_rearrangement.tsv` files. [ConvertDb fasta](https://changeo.readthedocs.io/en/stable/tools/ConvertDb.html#convertdb-py-fasta) is used to
generate a `.fasta` file from the rearrangement table.

### Assign genes with Igblast

<details markdown="1">
<summary>Output files</summary>

- `vdj_annotation/01-assigngenes/<sampleID>`
  - `*.fmt7`: Igblast results.
  - `*.fasta`: Igblast results converted to fasta.
  - `*log.txt`: Log of the process that will be parsed to generate a report.

</details>

Assign genes with Igblast, using the IMGT database is performed by the [AssignGenes](https://changeo.readthedocs.io/en/stable/examples/igblast.html#running-igblast) command of the Change-O tool from the Immcantation Framework.

### Make database from assigned genes

<details markdown="1">
<summary>Output files</summary>

- `vdj_annotation/02-makedb/<sampleID>`
  - `*log.txt`: Log of the process that will be parsed to generate a report.
  - `*db-pass.tsv`: Rearrangement table in AIRR-C format containing the assigned gene information.

</details>

IgBLAST's results are parsed and standardized with [MakeDB](https://changeo.readthedocs.io/en/stable/examples/igblast.html#processing-the-output-of-igblast) to follow the [AIRR Community standards](https://docs.airr-community.org/en/stable/datarep/rearrangements.html) for rearrangement data.

### Quality filter sequences

<details markdown="1">
<summary>Output files</summary>

- `vdj_annotation/03-quality-filter/<sampleID>`
  - `*log.txt`: Log of the process that will be parsed to generate a report.
  - `*quality-pass.tsv*`: Rearrangement table in AIRR-C format containing the sequences that passed the quality filtering steps.

</details>

A table is generated that retains sequences with concordant locus in the `v_call` and `locus` fields, with a `sequence_alignment` with a maximum of 10% of Ns and a length of at least 200 informative nucleotides (not `-`, `.` or `N`).

### Removal of non-productive sequences

<details markdown="1">
<summary>Output files</summary>

- `vdj_annotation/04-select-productive/<sampleID>`
  - `*log.txt`: Log of the process that will be parsed to generate a report.
  - `*productive-T.tsv*`: Rearrangement table in AIRR-C format, with only productive sequences.

</details>

Non-functional sequences identified with IgBLAST are removed with [ParseDb](https://changeo.readthedocs.io/en/stable/tools/ParseDb.html).

### Removal of sequences with junction length not multiple of 3

<details markdown="1">
<summary>Output files</summary>

- `vdj_annotation/05-select-junction-mod3/<sampleID>`
  - `*log.txt`: Log of the process that will be parsed to generate a report.
  - `*junction-pass.tsv*`: Rearrangement table in AIRR-C format, with only sequences that have a nucleotide junction length multiple of 3.

</details>

### Add metadata

<details markdown="1">
<summary>Output files</summary>

- `vdj_annotation/06-annotate-metadata/<sampleID>`
  - `*log.txt`: Log of the process that will be parsed to generate a report.
  - `*meta-pass.tsv*`: Rearrangement table in AIRR-C format annotated with metadata provided in the starting metadata sheet.

</details>

## Shazam

### Reconstruct germlines

<details markdown="1">
<summary>Output files</summary>

- `qc-filtering/bulk-qc-filtering/01-create-germlines/<sampleID>`
  - `*log.txt`: Log of the process that will be parsed to generate a report.
  - `*germ-pass.tsv`: Rearrangement table in AIRR-C format with an additional
    field with the reconstructed germline sequence for each sequence.

</details>

Reconstructing the germline sequences with the [CreateGermlines](https://changeo.readthedocs.io/en/stable/tools/CreateGermlines.html#creategermlines) Immcantation tool.

### Chimera filter

<details markdown="1">
<summary>Output files</summary>

- `qc-filtering/bulk-qc-filtering/02-chimera-filter/<sampleID>`
  - `*log.txt`: Log of the process that will be parsed to generate a report.
  - `*chimera-pass.tsv`: Rearrangement table in AIRR-C format sequences that
    passed the chimera removal filter.
  - `<sampleID>_chimera_report`: Report with plots showing the mutation patterns

</details>

Mutations patterns in different window sizes are analyzed with functions from
the Immcantation R package [SHazaM](https://shazam.readthedocs.io/en/stable/).

### Detect contamination

<details markdown="1">
<summary>Output files. Optional. </summary>

- `qc-filtering/bulk-qc-filtering/03-detect_contamination`
  - `*log.txt`: Log of the process that will be parsed to generate a report.
  - `*cont-flag.tsv`: Rearrangement table in AIRR-C format with sequences that
    passed the chimera removal filter.
  - `all_reps_cont_report`: Report.

</details>

This folder is genereated when `detect_contamination` is set to `true`.

### Collapse duplicates

<details markdown="1">
<summary>Output files. </summary>

- `qc-filtering/bulk-qc-filtering/04-collapse-duplicates/<sampleID>`
  - `*log.txt`: Log of the process that will be parsed to generate a report.
  - `*collapse_report/`: Report.
    - `repertoires/*collapse-pass.tsv`: Rearrangement table in AIRR-C format with duplicated
      sequences removed.

</details>

### Single cell QC

<details markdown="1">
<summary>Output files. </summary>

- `qc-filtering/single-cell-qc/all_reps_scqc_report`
  - `*log.txt`: Log of the process that will be parsed to generate a report.
  - `*all_reps_scqc_report/`: Report.
    - `*scqc-pass.tsv`: Rearrangement table in AIRR-C format with sequences that
      passed the quality filtering.

</details>

### Determining hamming distance threshold

<details markdown="1">
<summary>Output files</summary>

- `clonal_analysis/find-threshold/`
  - `*log`: Log of the process that will be parsed to generate a report.
  - `all_reps_threshold-mean.tsv`: Mean of all hamming distance thresholds of the
    Junction regions as determined by Shazam.
  - `all_reps_threshold-summary.tsv`: Thresholds for each group of `--cloneby` samples.
  - `all_reps_dist_report`: Report

</details>

Determining the hamming distance threshold of the junction regions for clonal determination using [Shazam](https://shazam.readthedocs.io) when `clonal_threshold` is set to `auto`.

## SCOPer define clones

### Define clones

<details markdown="1">
<summary>Output files</summary>

- `clonal_analysis/define_clones/<subjectID>`
  - `*log`: Log of the process that will be parsed to generate a report.
  - `repertoires/<sampleID>_clone-pass.tsv`: Rearrangement tables in AIRR-C format with sequences that
    passed the clonal assignment step. The field `clone_id` contains the clonal clusters identifiers.
  - `tables/`: Table in AIRR format containing the assigned gene information and an additional field with the clone id.
    - `clonal_abundance.tsv`
    - `clonal_diversity.tsv`
    - `clone_sizes_table.tsv`
    - `num_clones_table_nosingle.tsv`
    - `num_clones_table.tsv`
  - `ggplots/`: Diversity and abundance plots as `ggplot` objects.
  - `figures/`: Clone size, diversity and abundance `png` plots.

A similar output folder `clonal_analysis/define_clones/all_reps_clone_report` is generated for all data.

</details>

Assigning clones to the sequences obtained from IgBlast with the [scoper::hierarchicalClones](https://scoper.readthedocs.io/en/stable/topics/hierarchicalClones/) Immcantation tool.

#

## Lineage reconstruction

<details markdown="1">
<summary>Output files</summary>

- `clonal_analysis/dowser_lineages/`
  - `<sampleID>*log`: Log of the process that will be parsed to generate a report.
  - `<sample1ID>_dowser_report`: Report

</details>

Reconstructing clonal lineage with [IgPhyML](https://igphyml.readthedocs.io/en/stable/) and
[dowser](https://dowser.readthedocs.io/en/stable/topics/getTrees/) from the Immcantation toolset.

## Repertoire comparison

<details markdown="1">
<summary>Output files</summary>

- `repertoire_analysis/repertoire_comparison/`
  - `all_data.tsv`: AIRR format table containing the processed sequence information for all subjects.
  - `Abundance`: contains clonal abundance calculation plots and tables.
  - `Diversity`: contains diversity calculation plots and tables.
  - `V_family`: contains V gene and family distribution calculation plots and tables.
- `Airrflow_report.html`: Contains the repertoire comparison results in an html report form: Abundance, Diversity, V gene usage tables and plots. Comparison between treatments and subjects.

</details>

Calculation of several repertoire characteristics (diversity, abundance, V gene usage) for comparison between subjects, time points and cell populations. An Rmarkdown report is generated with the [Alakazam R package](https://alakazam.readthedocs.io/en/stable/).

## Tracking number of reads

<details markdown="1">
<summary>Output files</summary>

- `report_file_size/file_size_report`: Report summarizing the number of sequences after the most important pipeline steps.
  - `tables/*tsv`: Tables with the number of sequences at each processing step.

</details>

Parsing the logs from the previous processes. Summary of the number of sequences left after each of the most important pipeline steps.

## Log parsing

<details markdown="1">
<summary>Output files</summary>

- `parsed_logs/`
  - `sequences_table`: table summarizing of the number of sequences after the most important pipeline steps.

</details>

Parsing the logs from the previous processes. Summary of the number of sequences left after each of the most important pipeline steps.

## Databases

Copy of the downloaded IMGT database by the process `fetch_databases`, used for the gene assignment step.

## MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
