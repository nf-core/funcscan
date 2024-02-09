# nf-core/funcscan: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.5dev - [date]

### `Added`

- [#322](https://github.com/nf-core/funcscan/pull/322) Updated all modules: introduce environment.yml files. (by @jasmezz)
- [#324](https://github.com/nf-core/funcscan/pull/324) Removed separate DeepARG test profile because database download is now stable. (by @jasmezz)
- [#332](https://github.com/nf-core/funcscan/pull/332) & [#327](https://github.com/nf-core/funcscan/pull/327) Merged pipeline template of nf-core/tools version 2.12.1 (by @jfy133, @jasmezz)

### `Fixed`

### `Dependencies`

| Tool    | Previous version | New version |
| ------- | ---------------- | ----------- |
| DeepARG | 1.0.2            | 1.0.4       |
| DeepBGC | 0.1.30           | 0.1.31      |
| MultiQC | 1.15             | 1.18        |

### `Deprecated`

## v1.1.4 - [2023-11-07]

### `Added`

### `Fixed`

- [#306](https://github.com/nf-core/funcscan/pull/306) Added new parameter `annotation_prokka_retaincontigheaders` to allow prokka to retain the original contig headers/locus tag. (by @darcy220606)
- [#307](https://github.com/nf-core/funcscan/pull/307) Fixed stability of deepARG tests by using Zenodo copy of database (❤️ to Gustavo Arango and Liqing Zhang for uploading, fix by @jfy133)
- [#310](https://github.com/nf-core/funcscan/pull/310) Fixed error when supplying uncompressed input; added "fas" file extension for FASTA files. (by @tavareshugo)
- [#311](https://github.com/nf-core/funcscan/pull/311) Merged pipeline template of nf-core/tools version 2.10. (by @jasmezz)

### `Dependencies`

| Tool          | Previous version | New version |
| ------------- | ---------------- | ----------- |
| AMRFinderPlus | 3.10.42          | 3.11.18     |
| Bakta         | 1.7.0            | 1.8.2       |
| MultiQC       | 1.14             | 1.15        |

### `Deprecated`

- FastQC

## v1.1.3 - [2023-08-11]

### `Added`

- [#290](https://github.com/nf-core/funcscan/pull/290) Merged pipeline template of nf-core/tools version 2.9, updated references. (by @jfy133)
- [#285](https://github.com/nf-core/funcscan/pull/285) Use nf-validation for samplesheet checking and added support for `fna.gz` input FASTA files. (by @louperelo, @mirpedrol, @jfy133)
- [#295](https://github.com/nf-core/funcscan/pull/295) Add Prokka to MultiQC output. (by @louperelo)

### `Fixed`

- [#296](https://github.com/nf-core/funcscan/pull/296) Fixed empty output when saving prodigal annotations. (reported by @louperelo, fix by @jasmezz)
- [#297](https://github.com/nf-core/funcscan/pull/297) Added check for empty annotation files prior going into screening. (❤️ to @alexhbnr for requesting, added by @jfy133)
- [#299](https://github.com/nf-core/funcscan/pull/299) Fixed pigz error with symlinks in Pyrodigal. (by @jasmezz)
- [#300](https://github.com/nf-core/funcscan/pull/300) Fixed wrong Pyrodigal channels being submitted to antiSMASH. (reported by Till Bayer, fix by @jasmezz)
- [#302](https://github.com/nf-core/funcscan/pull/302) Removed trouble-causing default parameters in json schema. (by @robsyme)

### `Dependencies`

| Tool   | Previous version | New version |
| ------ | ---------------- | ----------- |
| comBGC | 0.6.0            | 0.6.1       |
| GECCO  | 0.9.2            | 0.9.8       |

### `Deprecated`

## v1.1.2 - [2023-06-30]

### `Added`

### `Fixed`

- [#279](https://github.com/nf-core/funcscan/pull/279) Fix docker/podman registry definition for tower compatibility. (♥️ to sunitj for reporting, fix by @adamrtalbot)

### `Dependencies`

### `Deprecated`

## v1.1.1 - [2023-05-24]

### `Added`

- [#270](https://github.com/nf-core/funcscan/pull/270) Merged pipeline template of nf-core/tools version 2.8 and updated modules accordingly. (by @jasmezz, @jfy133)
- [#274](https://github.com/nf-core/funcscan/pull/274) Update all modules: changed docker links according to the change of quay.io as default repository and Pyrodigal annotation output now zipped. (by @jasmezz)
- [#275](https://github.com/nf-core/funcscan/pull/275) Save DRAMP database in the common database directory if `--save_databases` is supplied. (by @jasmezz)

### `Fixed`

- [#272](https://github.com/nf-core/funcscan/pull/272) Fix typo in Prokka output path in modules.config. (by @jasmezz)
- [#273](https://github.com/nf-core/funcscan/pull/273) Update Ampir module after input bugfix in module. (reported by @mathavanpu, fix by @louperelo)
- [#276](https://github.com/nf-core/funcscan/pull/276) Fix Pyrodigal parameters in modules.config. (by @jasmezz)

### `Dependencies`

### `Deprecated`

## v1.1.0 - British Beans on Toast - [2023-04-27]

### `Added`

- [#238](https://github.com/nf-core/funcscan/pull/238) Added dedicated DRAMP database downloading step for AMPcombi to prevent parallel downloads when no database provided by user. (by @jfy133)
- [#235](https://github.com/nf-core/funcscan/pull/235) Added parameter `annotation_bakta_db_downloadtype` to be able to switch between downloading either full (33.1 GB) or light (1.3 GB excluding UPS, IPS, PSC, see parameter description) versions of the Bakta database. (by @jasmezz)
- [#249](https://github.com/nf-core/funcscan/pull/249) Added bakta annotation to CI tests. (by @jasmezz)
- [#251](https://github.com/nf-core/funcscan/pull/251) Added annotation tool: Pyrodigal. (by @jasmezz)
- [#252](https://github.com/nf-core/funcscan/pull/252) Added a new parameter `-arg_rgi_savejson` that saves the file `<samplename>.json` in the RGI directory. The default ouput for RGI is now only `<samplename>.txt`. (by @darcy220606)
- [#253](https://github.com/nf-core/funcscan/pull/253) Updated Prodigal to have compressed output files. (by @jasmezz)
- [#262](https://github.com/nf-core/funcscan/pull/262) Added comBGC function to screen whole directory of antiSMASH output (one subfolder per sample). (by @jasmezz)
- [#263](https://github.com/nf-core/funcscan/pull/263) Removed `AMPlify` from test_full.config. (by @jasmezz)
- [#266](https://github.com/nf-core/funcscan/pull/266) Updated README.md with Pyrodigal. (by @jasmezz)

### `Fixed`

- [#243](https://github.com/nf-core/funcscan/pull/243) Compress the ampcombi_complete_summary.csv in the output directory. (by @louperelo)
- [#237](https://github.com/nf-core/funcscan/pull/237) Reactivate DeepARG automatic database downloading and CI tests as server is now back up. (by @jfy133)
- [#235](https://github.com/nf-core/funcscan/pull/235) Improved annotation speed by switching off pipeline-irrelevant Bakta annotation steps by default. (by @jasmezz)
- [#235](https://github.com/nf-core/funcscan/pull/235) Renamed parameter `annotation_bakta_db` to `annotation_bakta_db_localpath`. (by @jasmezz)
- [#242](https://github.com/nf-core/funcscan/pull/242) Fixed MACREL '.faa' issue that was generated when it was run on its own and upgraded MACREL from version `1.1.0` to `1.2.0` (by @Darcy220606)
- [#248](https://github.com/nf-core/funcscan/pull/248) Applied best-practice `error("message")` to all (sub)workflow files. (by @jasmezz)
- [#254](https://github.com/nf-core/funcscan/pull/254) Further resource optimisation based on feedback from 'real world' datasets. (ongoing, reported by @alexhbnr and @Darcy220606, fix by @jfy133)
- [#266](https://github.com/nf-core/funcscan/pull/266) Fixed wrong process name in base.config. (reported by @Darcy220606, fix by @jasmezz)

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
