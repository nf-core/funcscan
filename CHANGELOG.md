# nf-core/funcscan: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.1.0 - [unreleased]

### `Added`

- [#421](https://github.com/nf-core/funcscan/pull/421) Updated to nf-core template 3.0.2 (by @jfy133)

### `Fixed`

### `Dependencies`

| Tool     | Previous version | New version |
| -------- | ---------------- | ----------- |
| AMPcombi | 0.2.2            | 2.0.1       |
| Macrel   | 1.2.0            | 1.4.0       |
| MultiQC  | 1.24.0           | 1.25.1      |

### `Deprecated`

## v2.0.0 - [2024-09-05]

### `Breaking change`

- [#391](https://github.com/nf-core/funcscan/pull/391) Made all "database" parameter names consistent, skip hmmsearch by default. (by @jasmezz)

| Old parameter                                    | New parameter                           |
| ------------------------------------------------ | --------------------------------------- |
| `annotation_bakta_db_localpath`                  | `annotation_bakta_db`                   |
| `arg_abricate_db`                                | `arg_abricate_db_id`                    |
| `arg_abricate_localdbdir`                        | `arg_abricate_db`                       |
| `arg_deeparg_data`                               | `arg_deeparg_db`                        |
| `arg_deeparg_data_version`                       | `arg_deeparg_db_version`                |
| `arg_rgi_database`                               | `arg_rgi_db`                            |
| `bgc_antismash_databases`                        | `bgc_antismash_db`                      |
| `bgc_antismash_installationdirectory`            | `bgc_antismash_installdir`              |
| `bgc_deepbgc_database`                           | `bgc_deepbgc_db`                        |
| `save_databases`                                 | `save_db`                               |
| `taxa_classification_mmseqs_databases_localpath` | `taxa_classification_mmseqs_db`         |
| `taxa_classification_mmseqs_databases_id`        | `taxa_classification_mmseqs_db_id`      |
| `taxa_classification_mmseqs_databases_savetmp`   | `taxa_classification_mmseqs_db_savetmp` |
| `amp_skip_hmmsearch`                             | `amp_run_hmmsearch`                     |
| `bgc_skip_hmmsearch`                             | `bgc_run_hmmsearch`                     |

- [#343](https://github.com/nf-core/funcscan/pull/343) Standardized the resulting workflow summary tables to always start with 'sample_id\tcontig_id\t..'. Reformatted the output of `hamronization/summarize` module. (by @darcy220606)
- [#411](https://github.com/nf-core/funcscan/pull/411) Optimised hAMRonization input: only high-quality hits from fARGene output are reported. (by @jasmezz, @jfy133)

### `Added`

- [#322](https://github.com/nf-core/funcscan/pull/322) Updated all modules: introduce environment.yml files. (by @jasmezz)
- [#324](https://github.com/nf-core/funcscan/pull/324) Removed separate DeepARG test profile because database download is now stable. (by @jasmezz)
- [#332](https://github.com/nf-core/funcscan/pull/332) & [#327](https://github.com/nf-core/funcscan/pull/327) Merged pipeline template of nf-core/tools version 2.12.1 (by @jfy133, @jasmezz)
- [#338](https://github.com/nf-core/funcscan/pull/338) Set `--meta` parameter to default for Bakta, with singlemode optional. (by @jasmezz)
- [#343](https://github.com/nf-core/funcscan/pull/343) Added contig taxonomic classification using [MMseqs2](https://github.com/soedinglab/MMseqs2/). (by @darcy220606)
- [#358](https://github.com/nf-core/funcscan/pull/358) Improved RGI databases handling, users can supply their own CARD now. (by @jasmezz)
- [#375](https://github.com/nf-core/funcscan/pull/375) Merged pipeline template of nf-core/tools version 2.14.1. (by @jfy133)
- [#381](https://github.com/nf-core/funcscan/pull/381) Added support for supplying pre-annotated sequences to the pipeline. (by @jfy133, @jasmezz)
- [#382](https://github.com/nf-core/funcscan/pull/382) Optimised BGC screening run time and prevent crashes due to too-short contigs by adding contig length filtering for BGC workflow only. (by @jfy133, @darcy220606)
- [#366](https://github.com/nf-core/funcscan/pull/366) Added nf-test on pipeline level. (by @jfy133, @Darcy220606, @jasmezz)
- [#403](https://github.com/nf-core/funcscan/pull/403) Added antiSMASH parameters `--pfam2go`, `--rre`, and `--tfbs`. (reported by @Darcy220606, added by @jasmezz)
- [#405](https://github.com/nf-core/funcscan/pull/405) Added argNorm to ARG subworkflow. (by @Vedanth-Ramji)

### `Fixed`

- [#348](https://github.com/nf-core/funcscan/pull/348) Updated samplesheet for pipeline tests to 'samplesheet_reduced.csv' with smaller datasets to reduce resource consumption. Updated prodigal module to fix pigz issue. Removed `tests/` from `.gitignore`. (by @darcy220606)
- [#362](https://github.com/nf-core/funcscan/pull/362) Save annotations from bakta in subdirectories per sample. (by @jasmezz)
- [#363](https://github.com/nf-core/funcscan/pull/363) Removed warning from DeepBGC usage docs. (by @jasmezz)
- [#365](https://github.com/nf-core/funcscan/pull/365) Fixed AMRFinderPlus module and usage docs for manual database download. (by @jasmezz)
- [#371](https://github.com/nf-core/funcscan/pull/371) Fixed AMRFinderPlus parameter `arg_amrfinderplus_name`. (by @m3hdad)
- [#377](https://github.com/nf-core/funcscan/pull/377) Fixed an occasional RGI process failure when certain files not produced. (❤️ to @amizeranschi for reporting, fix by @amizeranschi & @jfy133)
- [#386](https://github.com/nf-core/funcscan/pull/386) Updated DeepBGC module to fix output file names, separate annotation step for all BGC tools, add warning if no BGCs found, fix MultiQC reporting of annotation workflow. (by @jfy133, @jasmezz)
- [#393](https://github.com/nf-core/funcscan/pull/393) & [#397](https://github.com/nf-core/funcscan/pull/397) Fixed a docker/singularity only error appearing when running with conda. (❤️ to @ewissel for reporting, fix by @jfy33 & @jasmezz)
- [#391](https://github.com/nf-core/funcscan/pull/391) Skip hmmmsearch by default to not crash pipeline if user provides no HMM files, updated docs. (by @jasmezz)
- [#397](https://github.com/nf-core/funcscan/pull/397) Removed deprecated AMPcombi module, fixed variable name in BGC workflow, updated minor parts in docs (usage, parameter schema). (by @jasmezz)
- [#402](https://github.com/nf-core/funcscan/pull/402) Fixed BGC length calculation for antiSMASH hits by comBGC. (by @jasmezz)
- [#406](https://github.com/nf-core/funcscan/pull/406) Fixed prediction tools not being executed if annotation workflow skipped. (by @jasmezz)
- [#407](https://github.com/nf-core/funcscan/pull/407) Fixed comBGC bug when parsing multiple antiSMASH files. (by @jasmezz)
- [#409](https://github.com/nf-core/funcscan/pull/409) Fixed argNorm overwriting its output for DeepARG. (by @jasmezz, @jfy133)
- [#412](https://github.com/nf-core/funcscan/pull/412) Improve all pre-run database download documentation. (by @jfy133)

### `Dependencies`

| Tool          | Previous version | New version |
| ------------- | ---------------- | ----------- |
| AMPcombi      | 0.1.7            | 0.2.2       |
| AMPlify       | 1.1.0            | 2.0.0       |
| AMRFinderPlus | 3.11.18          | 3.12.8      |
| antiSMASH     | 6.1.1            | 7.1.0       |
| argNorm       | NA               | 0.5.0       |
| bioawk        | 1.0              | NA          |
| comBGC        | 1.6.1            | 1.6.2       |
| DeepARG       | 1.0.2            | 1.0.4       |
| DeepBGC       | 0.1.30           | 0.1.31      |
| GECCO         | 0.9.8            | 0.9.10      |
| hAMRonization | 1.1.1            | 1.1.4       |
| HMMER         | 3.3.2            | 3.4         |
| MMSeqs        | NA               | 2:15.6f452  |
| MultiQC       | 1.15             | 1.24        |
| Pyrodigal     | 2.1.0            | 3.3.0       |
| RGI           | 5.2.1            | 6.0.3       |
| seqkit        | NA               | 2.8.1       |
| tabix/htslib  | 1.11             | 1.20        |

### `Deprecated`

- [#384](https://github.com/nf-core/funcscan/pull/384) Deprecated AMPcombi and exchanged it with full suite of AMPcombi2 submodules. (by @darcy220606)
- [#382](https://github.com/nf-core/funcscan/pull/382) Optimised BGC screening run time and prevent crashes due to too-short contigs by adding contig length filtering for BGC workflow only. Bioawk is replaced with seqkit. (by @jfy133, @darcy220606)

## v1.1.6 - [2024-07-08]

### `Added`

### `Fixed`

- [#396](https://github.com/nf-core/funcscan/pull/396) Fixed bioawk overwriting input files. (❤️ to @Microbion for reporting, fix by @jfy133)

### `Dependencies`

## v1.1.5 - [2024-03-20]

### `Added`

### `Fixed`

- [#346](https://github.com/nf-core/funcscan/pull/346) Pinned version of nf-validation to 1.1.3

### `Dependencies`

| Plugin        | Previous | New version |
| ------------- | -------- | ----------- |
| Bakta         | 1.8.2    | 1.9.3       |
| nf-validation | Latest   | 1.1.3       |

### `Deprecated`

## v1.1.4 - [2023-11-07]

### `Added`

### `Fixed`

- [#306](https://github.com/nf-core/funcscan/pull/306) Added new parameter `annotation_prokka_retaincontigheaders` to allow prokka to retain the original contig headers/locus tag. (by @darcy220606)
- [#307](https://github.com/nf-core/funcscan/pull/307) Fixed stability of deepARG tests by using Zenodo copy of database. (❤️ to Gustavo Arango and Liqing Zhang for uploading, fix by @jfy133)
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
