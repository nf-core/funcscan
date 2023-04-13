# nf-core/funcscan: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - [unreleased]

### `Added`

- [#238](https://github.com/nf-core/funcscan/pull/238) Added dedicated DRAMP database downloading step for AMPcombi to prevent parallel downloads when no database provided by user (by @jfy133)
- [#235](https://github.com/nf-core/funcscan/pull/235) Added parameter `annotation_bakta_db_downloadtype` to be able to switch between downloading either full (33.1 GB) or light (1.3 GB excluding UPS, IPS, PSC, see parameter description) versions of the Bakta database. (by @jasmezz)
- [#252](https://github.com/nf-core/funcscan/pull/252) Added a new parameter `-arg_rgi_savejson` that saves the file `samplename.json` in the RGI directory. The default ouput for RGI is now only `samplename.txt`. (by @darcy220606)

### `Fixed`

- [#243](https://github.com/nf-core/funcscan/pull/243) Compress the ampcombi_complete_summary.csv in the output directory (by @louperelo)
- [#237](https://github.com/nf-core/funcscan/pull/237) Reactivate DeepARG automatic database downloading and CI tests as server is now back up. (by @jfy133)
- [#235](https://github.com/nf-core/funcscan/pull/235) Improved annotation speed by switching off pipeline-irrelevant Bakta annotation steps by default. (by @jasmezz)
- [#235](https://github.com/nf-core/funcscan/pull/235) Renamed parameter `annotation_bakta_db` to `annotation_bakta_db_localpath`. (by @jasmezz)
- [#248](https://github.com/nf-core/funcscan/pull/248) Applied best-practice `error("message")` to all (sub)workflow files. (by @jasmezz)

### `Dependencies`

| Tool  | Previous version | New version |
| ----- | ---------------- | ----------- |
| Bakta | 1.6.1            | 1.7.0       |

### `Deprecated`

## v1.0.1 - [2023-02-27]

### `Added`

- [#229](https://github.com/nf-core/funcscan/pull/229) Added pipeline DOI to `WorkflowMain.groovy` to display citation info when executing the pipeline. (by @jasmezz)

### `Fixed`

- [#227](https://github.com/nf-core/funcscan/pull/227) Removed a header check in the `check_samplesheet.py` script that was producing false negatives. Presence of required columns is still validated. (by @Midnighter)
- [#228](https://github.com/nf-core/funcscan/pull/228) Improved database downloading guidance to emphasise it is recommended to let nf-core/funcscan do the downloading on a first run, rather than manually downloading yourself. (reported by @alexhbnr, fixed by @jfy133)

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
