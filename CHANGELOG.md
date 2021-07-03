# nf-core/bcellmagic: Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [2.0.0dev] - date

* The pipeline was ported to the Nextflow DSL2 syntax.

### `Added`

* [#69](https://github.com/nf-core/bcellmagic/pull/69): Template update to nf-core tools v1.10.2
* [#69](https://github.com/nf-core/bcellmagic/pull/69): Added parameter json schema
* [#74](https://github.com/nf-core/bcellmagic/pull/74): Added possibility of setting UMI start position and better UMI docs
* [#85](https://github.com/nf-core/bcellmagic/pull/85): Added primer handling params: `--vprimer_start`, `--cprimer_start`, `--primer_mask_mode`, `--primer_maxeror`, `--primer_consensus`
* Template update to nf-core tools v1.12.1
* [#85](https://github.com/nf-core/bcellmagic/pull/85): Added support for TCR data with params: `--loci`
* [#85](https://github.com/nf-core/bcellmagic/pull/85): Added support for mice data and possibility for other species with params: `--species`
* [#85](https://github.com/nf-core/bcellmagic/pull/85): Added support for 5' RACE technology with params: `--race_5prime`, `--race_linker`
* [#87](https://github.com/nf-core/bcellmagic/pull/87): Added tests for TCR data
* [#102](https://github.com/nf-core/bcellmagic/pull/102): Added param `--protocol`.
* [#102](https://github.com/nf-core/bcellmagic/pull/102): Added support for C-primer in any R1 or R2 with param `--cprimer_position`.
* [#102](https://github.com/nf-core/bcellmagic/pull/102): Bump versions to 2.0 and NXF 21.04.0.
* [#103](https://github.com/nf-core/bcellmagic/pull/103): Added full size tests (pcr_umi).
* Added parameter `--skip_lineages`.
* [#112](https://github.com/nf-core/bcellmagic/pull/112): Template update to nf-core tools v1.14

### `Fixed`

* [#74](https://github.com/nf-core/bcellmagic/pull/74): Fixed AWStest workflow
* [#75](https://github.com/nf-core/bcellmagic/pull/75): Fixed lineage trees publish plots and graphml
* [#75](https://github.com/nf-core/bcellmagic/pull/75): Assemble process shorten filename to avoid error due to too long filename
* [#87](https://github.com/nf-core/bcellmagic/pull/87): Order sequence logs by sample ID.
* [#104](https://github.com/nf-core/bcellmagic/pull/104): Fix bug in pairseq barcode copy before consensus.
* Analysis not restricted to Ig heavy chains.

### `Dependencies`

### `Deprecated`

* [#69](https://github.com/nf-core/bcellmagic/pull/69): `--SkipDownstream` param changed to `--skip_downstream`.
* [#69](https://github.com/nf-core/bcellmagic/pull/69): `--metadata` param changed to `--input`
* `--saveDBs` param change to `--save_databases`
* Default for `--umi_length` changed to 0
* [#102](https://github.com/nf-core/bcellmagic/pull/102): `--race_5prime` param deprecated in favor of `--protocol`.
* `--skip_downstream` changed to `--skip_report`.

## [1.2.0] - 2020-01-14 - "Riddikulus"

### `Added`

* Handle barcodes that are already merged to R1 or R2 reads
* Validate inputs and cluster threshold
* `--downstream_only` feature
* Handle of UMIs of different lengths
* `--skipDownstream` feature
* Add github actions ci testing

### `Fixed`

* [#51](https://github.com/nf-core/bcellmagic/issues/51) - Fixed MaskPrimers bug
* [#45](https://github.com/nf-core/bcellmagic/issues/45) - Fixed UMI reading from R1 or R2 & UMI length
* [#57](https://github.com/nf-core/bcellmagic/issues/57) - Improved results directory organization
* [#55](https://github.com/nf-core/bcellmagic/issues/55) - Dropped Singularity file

### `Dependencies`

### `Deprecated`

## [1.1.0] - 2019-11-06 - "Wingardium Leviosa"

### `Added`

* Merging all the repertoires from the same patient
* Added clone calculation per patient: `--set_cluster_threshold`and `cluster_threshold` parameters.
* Added downstream analysis processes: diversity, abundance, mutational load, Ig type and gene distribution
* Parsing logs for all processes
* FastQC and multiQC processes
* Option for providing Illumina index and UMI as part of R1
* Update template to tools `1.7`

### `Fixed`

* [#25](https://github.com/nf-core/bcellmagic/issues/25) - Improved documentation
* [#27](https://github.com/nf-core/bcellmagic/issues/27) - Added FastQC and MultiQC processes
* [#21](https://github.com/nf-core/bcellmagic/issues/21) - Added log parsing

### `Dependencies`

* Update Nextflow `0.32.0` -> `19.10.0`
* Added several requirements for downstream analysis.

### `Deprecated`

## [1.0.0] - 2019-04-16 - "Alohomora"

* Initial release of nf-core/bcellmagic, created with the [nf-core](http://nf-co.re/) template.
