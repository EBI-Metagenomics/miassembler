# ebi-metagenomics/miassembler: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v3.1.15 - 2026-09-01

### Fixed

- Fixed a typo that affected metaQuast when processing metaT data ([#85](https://github.com/EBI-Metagenomics/miassembler/pull/85))

## v3.1.14 - 2026-06-11

### Changed

- Upgrade bwa-mem2 to version 2.3 ([#84](https://github.com/EBI-Metagenomics/miassembler/pull/84))

## v3.1.13 - 2026-04-08

### Fixed

- Missing join before `minimap_align` in LR / LQ subworkflows ([#82](https://github.com/EBI-Metagenomics/miassembler/pull/82))
- Add join for frameshift correction input ([#82](https://github.com/EBI-Metagenomics/miassembler/pull/82))

### Changed

- Refactor time and memory resource allocation for long reads ([#81](https://github.com/EBI-Metagenomics/miassembler/pull/81))
- Remove commented-out process for publishing cleaned contigs ([#81](https://github.com/EBI-Metagenomics/miassembler/pull/81))

## v1.0dev - [date]

Initial release of ebi-metagenomics/miassembler, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
