# nf-core/funcscan: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev [unreleased]

### `Added`

- Added pipeline citation DOI to `WorkflowMain.groovy` to display it when executing the pipeline

### `Fixed`

- Removed a header check in the `check_samplesheet.py` script that was producing false negatives. Presence of required columns is still validated. (by @Midnighter)
- Improved database downloading guidance to emphasise it is recommended to let nf-core/funcscan do the downloading on a first run, rather than manually downloading yourself (reported by @alexhbnr, fixed by @jfy133)

### `Dependencies`

### `Deprecated`

## v1.0.0 - German Rollmops - [2023-02-15]

Initial release of nf-core/funcscan, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Added annotation tools (Prokka, Prodigal, Bakta).
- Added AMP screening workflow (tools: Macrel, AMPlify, ampir, hmmsearch).
- Added ARG screening workflow (tools: ABRicate, AMRFinderPlus, DeepARG, fARGene).
- Added BGC screening workflow (tools: antiSMASH, DeepBGC, GECCO, hmmsearch).
- Added workflow summary tools (tools: hAMRonization, AMPcombi, comBGC).

### `Fixed`

### `Dependencies`

### `Deprecated`
