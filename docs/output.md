# nf-core/bcellmagic: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

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
  * Fastq with only reads that passed the quality filter.
* `fastq/*.tab`
  * table containing read ID and quality.
* `fastq/*.log`
  * Log of the process.

## Mask primers
Masks primers that are provided in the C-primers and V-primers input files. It uses the tool [MaskPrimers](https://presto.readthedocs.io/en/version-0.5.11/tools/MaskPrimers.html) of the Presto Immcantation toolset.

**Output directory: `results/mask_primers`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `fastq/*.fastq`
  * Fastq with reads with masked primers.
* `fastq/*.log`
  * Log containing sequence identifiers and the error in masking primers.

## Pair sequences
Pair read mates using [PairSeq](https://presto.readthedocs.io/en/version-0.5.11/tools/PairSeq.html) from the Presto Immcantation toolset.

**Output directory: `results/pair_sequences`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `fastq/*.fastq`
  * Fastq with reads that passed mate pairing.

## Cluster sets
Cluster sequences according to similarity, using [ClusterSets set](https://presto.readthedocs.io/en/version-0.5.11/tools/ClusterSets.html#clustersets-set). This step is introduced to deal with too low UMI diversity.

**Output directory: `results/cluster_sets`**

* `command_log.txt``
  * Log of the process that will be parsed to generate a report.
* `fastq/*.fastq`
  * Fastq with reads and annotation in their headers of cluster group.

## Build UMI consensus
Build consensus of UMI from all sequences that were annotated to have the same UMI. Uses [BuildConsensus](https://presto.readthedocs.io/en/version-0.5.11/tools/BuildConsensus.html).

**Output directory: `results/build_consensus`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `fastq/*.fastq`
  * Fastq with reads that passed the build consensus ste.
* `info/*.tab`
  * Parsed log containing the sequence barcodes and primers info

## Re-pair mates
Re-pair read mates using [PairSeq](https://presto.readthedocs.io/en/version-0.5.11/tools/PairSeq.html) from the Presto Immcantation toolset.

**Output directory: `results/repair_mates`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `fastq/*.fastq`
  * Fastq with reads that passed mate pairing.
* `info/*.tab`
  * Parsed log contaning the sequence barcodes and re-pair info.

## Assemble mates
Assemble read mates using [AssemblePairs](https://presto.readthedocs.io/en/version-0.5.11/tools/AssemblePairs.html) from the Presto Immcantation toolset.

**Output directory: `results/assemble_pairs`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `fastq/*.fastq`
  * Fastq with assembled reads.
* `info/*.tab`
  * Parsed log contaning the sequence barcodes and assemble pairs.

## Remove duplicates
Remove duplicates using [CollapseSeq](https://presto.readthedocs.io/en/version-0.5.11/tools/CollapseSeq.html) from the Presto Immcantation toolset.

**Output directory: `results/deduplicates`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `fastq/*.fastq`
  * Fastq with de-duplicated reads.
* `info/*.tab`
  * Parsed log contaning the sequence barcodes and deduplicated pairs.

## Filter sequences for at least 2 representative
Remove duplicates using [SplitSeq](https://presto.readthedocs.io/en/version-0.5.11/tools/SplitSeq.html) from the Presto Immcantation toolset.

**Output directory: `results/filter_representative_2`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `fastq/*.fastq`
  * Fastq with reads that have at least 2 representatives.
* `info/*.tab`
  * Parsed log contaning the sequence barcodes and split seq information.

## Assign genes with IG Blast
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
* `table/igh_genotyped_clone-pass.tab`
  * Table in ChangeO format contaning the assigned gene information and an additional field with the clone number.
* `info/igh_genotyped_table.tab`
  * Parsed log with sequence ID, assigned gene calls, junction length and clones.

## Reconstructing germlines
Reconstructing the germline sequences with the [CreateGermlines](https://changeo.readthedocs.io/en/version-0.4.5/tools/CreateGermlines.html#creategermlines) Immcantation tool.

**Output directory: `results/define_clones`**

* `command_log.txt`
  * Log of the process that will be parsed to generate a report.
* `table/igh_genotyped_clone-pass_germ-pass.tab`
  * Table in ChangeO format contaning the assigned gene information and an additional field with the germline reconstructed gene calls.

## Alakazam - repertoire analysis
Repertoire analysis with the Alakazam R package from the Immcantation toolset.

**Output directory: `results/alakazam`**

* `igh_genotyped_clone-pass_germ-pass.tab`
  * Final table in ChangeO format contaning the assigned gene information and an additional field with the germline reconstructed gene calls.



