# nf-core/bcellmagic: Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [1.3.0dev] - date

### `Added`

* [#69](https://github.com/nf-core/bcellmagic/pull/69): Template update to nf-core tools v1.10.2
* [#69](https://github.com/nf-core/bcellmagic/pull/69): Added parameter json schema
* [#74](https://github.com/nf-core/bcellmagic/pull/74): Added possibility of setting UMI start position and better UMI docs
* Added primer handling params: `--vprimer_start`, `--cprimer_start`, `--primer_mask_mode`, `--primer_maxeror`, `--primer_consensus`
* Template update to nf-core tools v1.12.1
* Added support for TCR data with params: `--loci`
* Added support for mice data and possibility for other species with params: `--species`
* Added support for 5' RACE technologiy with params: `--race_5prime`, `--race_linker`

### `Fixed`

* [#74](https://github.com/nf-core/bcellmagic/pull/74): Fixed AWStest workflow
* [#75](https://github.com/nf-core/bcellmagic/pull/75): Fixed lineage trees publish plots and graphml
* [#75](https://github.com/nf-core/bcellmagic/pull/75): Assemble process shorten filename to avoid error due to too long filename

### `Dependencies`

### `Deprecated`

* [#69](https://github.com/nf-core/bcellmagic/pull/69): `--SkipDownstream` param changed to `--skip_downstream`
* [#69](https://github.com/nf-core/bcellmagic/pull/69): `--metadata` param changed to `--input`
* `--saveDBs` param change to `--save_databases`
* Default for `--umi_length` changed to 0

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
