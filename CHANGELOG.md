# nf-core/airrflow: Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [2.2.1] - dev

### `Added`

- [#180](https://github.com/nf-core/airrflow/pull/180) Added possibility to add any property in the AIRR sequence table as label on the lineage tree nodes.
- [#180](https://github.com/nf-core/airrflow/pull/180) Lineage tree construction now also includes trees with just one sequence.
- [#180](https://github.com/nf-core/airrflow/pull/180) Added metadata annotation to final repertoire table.
- [#180](https://github.com/nf-core/airrflow/pull/180) Added possibility to provide custom Rmarkdown report.
- [#183](https://github.com/nf-core/airrflow/pull/183) Update template to nf-core tools v2.5.2

### `Fixed`

- [#180](https://github.com/nf-core/airrflow/pull/180) Repertoire analysis report now also saves diversity table.

### `Dependencies`

### `Deprecated`

## [2.2.0] - 2022-06-02 "Reparo"

### `Added`

- Pulling IMGT database cache from aws for CI tests.
- Added test to pull database from IMGT and build it with igblast.
- Template update to nf-core tools v2.4.1.
- Added zipped DB cache to nf-core/test-datasets.

### `Fixed`

- Updated container of `Fetch databases` and `Changeo_assigngenes` process to `Changeo=1.2.0` and `Igblast=1.17.1`and extended biocontainers base to have internet access.
- Fixed publishing directory mode for all modules.

### `Dependencies`

### `Deprecated`

## [2.1.0] - 2022-05-02 "Accio"

### `Added`

- [#130](https://github.com/nf-core/bcellmagic/pull/130): Organized presto processes in `presto_umi` subworkflow.
- [#128](https://github.com/nf-core/bcellmagic/pull/128): Added `presto_sans_umi` subworkflow option. Added postassembly FastQC and corresponding section in MultiQC. Included refs for analysis of light chains (if present) by default.
- updated docs for `--library_generation_method` parameter.
- Samplesheet column names were updated to follow the AIRR standard.
- Fixed docs on `--umi_start` parameter, this parameter should only be used when UMIs are provided in the index reads.
- [#143](https://github.com/nf-core/airrflow/pull/128) Template update to nf-core tools v2.2.
- [#150](https://github.com/nf-core/airrflow/pull/150) Added option to search for reverse primers.
- [#159](https://github.com/nf-core/airrflow/pull/159) Template update to nf-core tools v2.3.1, v2.3.1
- [#161](https://github.com/nf-core/airrflow/pull/161) Add option to skip clustering sequences in the UMI workflow
- [#163](https://github.com/nf-core/airrflow/pull/163) Added process labels and software version emitting to all modules. Fixed output folder name for changeo processes.
- Bump versions to 2.1.0

### `Fixed`

- [#150](https://github.com/nf-core/airrflow/pull/150): Fixed cprimer start position, when cprimer in R2 reads.
- Remove need for a plot when Hamming threshold cannot be generated.
- The shazam threshold process is not executed when the hamming threshold is provided.
- [#164](https://github.com/nf-core/airrflow/pull/150): Fixed AWS tests when running on fusion mounts, solving [#137](https://github.com/nf-core/airrflow/issues/137).

### `Dependencies`

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |

| igblast | 1.15.0 | 1.17.1 |
| presto | 0.6.2 | 0.7.0 |
| changeo | 1.0.2 | 1.2.0 |
| r-base | 4.0.3 | 4.1.2 |
| r-alakazam | 1.0.2 | 1.2.0 |
| r-shazam | 0.1.11 | 1.1.0 |
| r-tigger | 0.3.1 | 1.0.0 |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if new version information isn't present.

### `Deprecated`

- Parameter `--protocol` was updated to `--library_generation_method` to follow the AIRR standard.
- Parameter `--loci` and `--species` were converted to a column in the samplesheet (`pcr_target_locus` and `species`) to allow processing simultaneously TR and IG loci from the same sample, also allow processing different species in one samplesheet.
- Removed genotyping step with tiGGer.

## [2.0.0] - 2021-07-19 "Lumos"

### :warning: Major enhancements

- The pipeline has been ported to Nextflow DSL2
- Analysis of TCR repertoires is now also supported
- Added Rmarkdown report with a summary of the repertoire analysis
- Analysis of data with RACE 5' protocol is now supported
- Scripts to download and build new versions of the IMGT database have been added
- Improvement of UMI handling and possibility to exchange position of C and V-primers
- Updated to new versions of Immcantation Framework, providing support for the AIRR format

### `Added`

- [#69](https://github.com/nf-core/bcellmagic/pull/69): Template update to nf-core tools v1.10.2
- [#69](https://github.com/nf-core/bcellmagic/pull/69): Added parameter json schema
- [#74](https://github.com/nf-core/bcellmagic/pull/74): Added possibility of setting UMI start position and better UMI docs
- [#85](https://github.com/nf-core/bcellmagic/pull/85): Added primer handling params: `--vprimer_start`, `--cprimer_start`, `--primer_mask_mode`, `--primer_maxeror`, `--primer_consensus`
- [#85](https://github.com/nf-core/bcellmagic/pull/85): Added support for TCR data with params: `--loci`
- [#85](https://github.com/nf-core/bcellmagic/pull/85): Added support for mice data and possibility for other species with params: `--species`
- [#85](https://github.com/nf-core/bcellmagic/pull/85): Added support for 5' RACE technology with params: `--protocol`, `--race_linker`
- [#87](https://github.com/nf-core/bcellmagic/pull/87): Added tests for TCR data
- [#102](https://github.com/nf-core/bcellmagic/pull/102): Added param `--protocol`.
- [#102](https://github.com/nf-core/bcellmagic/pull/102): Added support for C-primer in any R1 or R2 with param `--cprimer_position`.
- [#102](https://github.com/nf-core/bcellmagic/pull/102): Bump versions to 2.0.
- [#103](https://github.com/nf-core/bcellmagic/pull/103): Added full size tests (pcr_umi).
- [#114](https://github.com/nf-core/bcellmagic/pull/114): Added parameter `--skip_lineages`.
- [#112](https://github.com/nf-core/bcellmagic/pull/112): Template update to nf-core tools v1.14.
- [#114](https://github.com/nf-core/bcellmagic/pull/114): Added Bcellmagic html report.
- [#114](https://github.com/nf-core/bcellmagic/pull/114): Improved documentation on amplicon protocol support.
- [#115](https://github.com/nf-core/bcellmagic/pull/115): Improved output file structure and documentation.
- [#124](https://github.com/nf-core/bcellmagic/pull/124): Template update to nf-core tools v2.0.1

### `Fixed`

- [#74](https://github.com/nf-core/bcellmagic/pull/74): Fixed AWStest workflow
- [#75](https://github.com/nf-core/bcellmagic/pull/75): Fixed lineage trees publish plots and graphml
- [#75](https://github.com/nf-core/bcellmagic/pull/75): Assemble process shorten filename to avoid error due to too long filename
- [#87](https://github.com/nf-core/bcellmagic/pull/87): Order sequence logs by sample ID.
- [#104](https://github.com/nf-core/bcellmagic/pull/104): Fix bug in pairseq barcode copy before consensus.
- [#114](https://github.com/nf-core/bcellmagic/pull/114): Analysis not restricted to Ig heavy chains.
- [#123](https://github.com/nf-core/bcellmagic/pull/123): Fix report Rmarkdown reading for running on AWS.

### `Dependencies`

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency   | Old version | New version |
| ------------ | ----------- | ----------- |
| Python       | v3.6.11     | 3.8.0,3.0.3 |
| markdown     | v3.1.1      |             |
| muscle       | v3.8.1551   |             |
| vsearch      | v2.11.1     |             |
| fastqc       | 0.11.8      | 0.11.9      |
| multiqc      | 1.9         | 1.11        |
| matplotlib   | 3.0.3       |             |
| cd-hit       | 4.8.1       |             |
| blast        | 2.7.1       |             |
| igblast      | 1.10.0      | 1.15.0      |
| phylip       | 3.697       | 3.697       |
| airr         | 1.2.1       |             |
| presto       | 0.5.10      | 0.6.2       |
| changeo      | 0.4.5       | 1.0.2       |
| biopython    | 1.70        | 1.74        |
| pandas       | 0.24.2      | 1.1.5       |
| seaborn      | 0.9.0       |             |
| r-base       | 3.5.1       | 4.0.3       |
| r-alakazam   | 0.2.11      | 1.0.2       |
| r-shazam     | 0.1.11      | 1.0.2       |
| r-tigger     | 0.3.1       | 1.0.0       |
| r-dplyr      | 0.8.3       | 1.0.6       |
| r-data.table | 1.12.2      |             |
| r-igraph     | 1.2.4.1     |             |
| igphyml      |             | 1.1.3       |
| wget         |             | 1.20.1      |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

### `Deprecated`

- [#69](https://github.com/nf-core/bcellmagic/pull/69): `--SkipDownstream` param changed to `--skip_downstream`.
- [#69](https://github.com/nf-core/bcellmagic/pull/69): `--metadata` param changed to `--input`
- `--saveDBs` param change to `--save_databases`
- Default for `--umi_length` changed to 0
- [#102](https://github.com/nf-core/bcellmagic/pull/102): `--race_5prime` param deprecated in favor of `--protocol`.
- `--skip_downstream` changed to `--skip_report`.

## [1.2.0] - 2020-01-14 - "Riddikulus"

### `Added`

- Handle barcodes that are already merged to R1 or R2 reads
- Validate inputs and cluster threshold
- `--downstream_only` feature
- Handle of UMIs of different lengths
- `--skipDownstream` feature
- Add github actions ci testing

### `Fixed`

- [#51](https://github.com/nf-core/bcellmagic/issues/51) - Fixed MaskPrimers bug
- [#45](https://github.com/nf-core/bcellmagic/issues/45) - Fixed UMI reading from R1 or R2 & UMI length
- [#57](https://github.com/nf-core/bcellmagic/issues/57) - Improved results directory organization
- [#55](https://github.com/nf-core/bcellmagic/issues/55) - Dropped Singularity file

### `Dependencies`

### `Deprecated`

## [1.1.0] - 2019-11-06 - "Wingardium Leviosa"

### `Added`

- Merging all the repertoires from the same patient
- Added clone calculation per patient: `--set_cluster_threshold`and `cluster_threshold` parameters.
- Added downstream analysis processes: diversity, abundance, mutational load, Ig type and gene distribution
- Parsing logs for all processes
- FastQC and multiQC processes
- Option for providing Illumina index and UMI as part of R1
- Update template to tools `1.7`

### `Fixed`

- [#25](https://github.com/nf-core/bcellmagic/issues/25) - Improved documentation
- [#27](https://github.com/nf-core/bcellmagic/issues/27) - Added FastQC and MultiQC processes
- [#21](https://github.com/nf-core/bcellmagic/issues/21) - Added log parsing

### `Dependencies`

- Update Nextflow `0.32.0` -> `19.10.0`
- Added several requirements for downstream analysis.

### `Deprecated`

## [1.0.0] - 2019-04-16 - "Alohomora"

- Initial release of nf-core/bcellmagic, created with the [nf-core](http://nf-co.re/) template.
