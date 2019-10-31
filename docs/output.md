# nf-core/bcellmagic: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [Filter sequence quality](#filter-sequence-quality) - filter sequences by quality
* [Mask primers](#mask-primers) - Masking primers
* [Pair mates](#pair-mates) - Pairing sequence mates.
* [Cluster sets](#cluster-sets) - Cluster sequences according to similarity.
* [Build consensus](#build-UMI-consensus) - Build UMI consensus
* [Re-pair mates](#re-pair-mates) - Re-pairing sequence mates.
* [Assemble mates](#assemble-mates) - Assemble sequence mates.
* [Remove duplicates](#remove-duplicates) - Remove read duplicates.
* [Filter sequences for at least 2 representative](#filter-sequences-for-at-least-2-representative) Filter sequences that do not have at least 2 reads assigned.
* [Assign genes with IgBlast](#assign-genes-with-igblast)
* [Determining genotype and hamming distance threshold](#determining-genotype-and-hamming-distance-threshold)
* [Defining clones](#defining-clones) - Defining clonal B-cell populations
* [Reconstructing germlines](#reconstructing-germlines) - Reconstruct gene calls of germline sequences
* [Clonal analysis](#clonal-analysis) - Clonal analysis.
* [Repertoire comparison](#repertoire-comparison) - Repertoire comparison.
* [Log parsing](#log-parsing) - Log parsing.
* [MultiQC](#MultiQC) - MultiQC

## FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

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
* `*.tab`
  * table containing read ID and quality.

## Mask primers
Masks primers that are provided in the C-primers and V-primers input files. It uses the tool [MaskPrimers](https://presto.readthedocs.io/en/version-0.5.11/tools/MaskPrimers.html) of the Presto Immcantation toolset.

**Output directory: `results/mask_primers`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `*.log`
  * Log containing sequence identifiers and the error in masking primers.

## Pair mates
Pair read mates using [PairSeq](https://presto.readthedocs.io/en/version-0.5.11/tools/PairSeq.html) from the Presto Immcantation toolset.

**Output directory: `results/pair_sequences`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.

## Cluster sets
Cluster sequences according to similarity, using [ClusterSets set](https://presto.readthedocs.io/en/version-0.5.11/tools/ClusterSets.html#clustersets-set). This step is introduced to deal with too low UMI diversity.

**Output directory: `results/cluster_sets`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.

## Build UMI consensus
Build consensus of UMI from all sequences that were annotated to have the same UMI. Uses [BuildConsensus](https://presto.readthedocs.io/en/version-0.5.11/tools/BuildConsensus.html).

**Output directory: `results/build_consensus`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `*.tab`
  * Parsed log containing the sequence barcodes and primers info

## Re-pair mates
Re-pair read mates using [PairSeq](https://presto.readthedocs.io/en/version-0.5.11/tools/PairSeq.html) from the Presto Immcantation toolset.

**Output directory: `results/repair_mates`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `*.tab`
  * Parsed log contaning the sequence barcodes and re-pair info.

## Assemble mates
Assemble read mates using [AssemblePairs](https://presto.readthedocs.io/en/version-0.5.11/tools/AssemblePairs.html) from the Presto Immcantation toolset.

**Output directory: `results/assemble_pairs`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `*.tab`
  * Parsed log contaning the sequence barcodes and assemble pairs.

## Remove duplicates
Remove duplicates using [CollapseSeq](https://presto.readthedocs.io/en/version-0.5.11/tools/CollapseSeq.html) from the Presto Immcantation toolset.

**Output directory: `results/deduplicates`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `*.tab`
  * Parsed log contaning the sequence barcodes and deduplicated pairs.

## Filter sequences for at least 2 representative
Remove sequences which do not have 2 representative using [SplitSeq](https://presto.readthedocs.io/en/version-0.5.11/tools/SplitSeq.html) from the Presto Immcantation toolset.

**Output directory: `results/filter_representative_2`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `*.tab`
  * Parsed log contaning the sequence barcodes and split seq information.

## Assign genes with IgBlast
Assign genes from the IGblast database using [AssignGenes](https://changeo.readthedocs.io/en/version-0.4.5/examples/igblast.html#running-igblast) and generating a table with [MakeDB](https://changeo.readthedocs.io/en/version-0.4.5/examples/igblast.html#processing-the-output-of-igblast). Non-functional sequences are removed with [ParseDb](https://changeo.readthedocs.io/en/version-0.4.5/tools/ParseDb.html). Sequences in are additionally converted to a fasta file with the [ConvertDb](https://changeo.readthedocs.io/en/version-0.4.5/tools/ConvertDb.html?highlight=convertdb) tool.

**Output directory: `results/igblast`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `fasta/*.fasta`
  * Blast results converted to fasta fall with genotype V-call annotated in the header.
* `table/*.tab`
  * Table in ChangeO format contaning the assigned gene information and metadata provided in the starting metadata sheet.

## Determining genotype and hamming distance threshold
Determining genotype and the hamming distance threshold of the junction regions for clonal determination using the [tigGER](https://tigger.readthedocs.io/en/0.3.1/) and [Shazam](https://shazam.readthedocs.io/en/version-0.1.11_a/).

**Output directory: `results/shazam`**

* `threshold.txt`
  * Hamming distance threshold of the Junction regions as determined by Shazam.
* `Hamming_distance_threshold.pdf`
  * Plot of the Hamming distance distribution between junction regions displaying the threshold for clonal assignment as determined by Shazam.
* `genotype.pdf`
  * Plot representing the patient genotype assessed by TigGER.
* `igh_genotyped.tab`
  * Table in ChangeO additionally containing the assigned genotype in V_CALL_GENOTYPED.
* `v_genotype.fasta`
  * Fasta file containing the full sequences for all V genes assigned to the patient.

## Defining clones
Assigning clones to the sequences obtained from IgBlast with the [DefineClones](https://changeo.readthedocs.io/en/version-0.4.5/tools/DefineClones.html?highlight=DefineClones) Immcantation tool.

**Output directory: `results/define_clones`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `igh_genotyped_clone-pass.tab`
  * Table in ChangeO format contaning the assigned gene information and an additional field with the clone number.
* `igh_genotyped_table.tab`
  * Parsed log with sequence ID, assigned gene calls, junction length and clones.

## Reconstructing germlines
Reconstructing the germline sequences with the [CreateGermlines](https://changeo.readthedocs.io/en/version-0.4.5/tools/CreateGermlines.html#creategermlines) Immcantation tool.

**Output directory: `results/define_clones`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `table/igh_genotyped_clone-pass_germ-pass.tab`
  * Table in ChangeO format contaning the assigned gene information and an additional field with the germline reconstructed gene calls.

## Clonal analysis
Reconstructing clonal linage with the Alakazam R package from the Immcantation toolset. Calculating and plotting several clone statistics.

**Output directory: `results/clonal_analysis`**

* `Clone_lineage/`
  * `Clones_table_patient.tsv`: contains a summary of the clones found for the patient, and the number of unique and total sequences identified in each clone.
  * `Clone_tree_plots`: contain a rooted graphical representation of each of the clones.
  * `Clone_lineage`: contain a GraphmL exported format of the plots. `All_graphs_patient.graphml` contains all graphs for that patient.
* `Clone_numbers/`
  * Number of clones and number of sequences per clone, patient-wise and cell population wise.
* `Clone_overlap/`
  * Plots for representing the clone overlap in number of clones and number of sequences between different time-points and cell populations of one patient.

## Repertoire comparison
Calculation of several repertoire characteristics (diversity, abundance) for comparison between patients, time points and cell popultions.

**Output directory: `results/repertoire_comparison`**

* `diversity/`
  * Diversity calculation
* `abundance/`
  * Abundance calculation
* `mutational_load/`
  * Mutational load

## Log parsing
Parsing the logs from the previous processes. 

**Output directory: `results/parsing_logs`**

* A table summarizing of the number of sequences after the most important steps is shown.

## MultiQC
[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)
