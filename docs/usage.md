# nf-core/funcscan: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/funcscan/usage](https://nf-co.re/funcscan/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

nf-core/funcscan is a pipeline for efficient and parallelised screening of long nucleotide sequences such as contigs for antimicrobial peptide genes, antimicrobial resistance genes, and biosynthetic gene clusters. It can additionally identify the taxonomic origin of the sequences.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/funcscan --input samplesheet.csv --outdir <OUTDIR> -profile docker --run_<amp/arg/bgc>_screening
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

To run any of the three screening workflows (AMP, ARG, and/or BGC) or taxonomic classification, switch them on by adding the respective flag(s) to the command:

- `--run_amp_screening`
- `--run_arg_screening`
- `--run_bgc_screening`
- `--run_taxa_classification`

When switched on, all tools of the given workflow will be run by default. If you don't need specific tools, you can explicitly skip them. For the taxonomic classification, MMseqs2 is currently the only tool implemented in the pipline.

**Example:** You want to run AMP and ARG screening but you don't need the DeepARG tool of the ARG workflow and the Macrel tool of the AMP workflow. Your command would be:

```bash
nextflow run nf-core/funcscan --input samplesheet.csv --outdir <OUTDIR> -profile docker --run_arg_screening --arg_skip_deeparg --run_amp_screening --amg_skip_macrel

```

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing temporary files required for the run
<OUTDIR>        # Final results (location specified with --outdir)
.nextflow_log   # Log file from nextflow

# Other nextflow hidden files, eg. history of pipeline runs and old logs
```

## Samplesheet input

nf-core/funcscan takes FASTA files as input, typically contigs or whole genome sequences. To supply these to the pipeline, you will need to create a samplesheet with information about the samples you would like to analyse. Use this parameter to specify its location.

```bash
--input '[path to samplesheet file]'
```

The input samplesheet has to be a comma-separated file (`.csv`) with 2 (`sample`, and `fasta`) or 4 columns (`sample`, `fasta`, `protein`, `gbk`), and a header row as shown in the examples below.

If you already have annotated contigs with peptide sequences and an annotation file in Genbank format (`.gbk.` or `.gbff`), you can supply these to the pipeline using the optional `protein` and `gbk` columns. If these additional columns are supplied, pipeline annotation (i.e. with bakta, prodigal, pyrodigal or prokka) will be skipped and the corresponding annotation files used instead.

For two columns (without pre-annotated data):

```csv title="samplesheet.csv"
sample,fasta
sample_1,/<path>/<to>/wastewater_metagenome_contigs_1.fasta.gz
sample_2,/<path>/<to>/wastewater_metagenome_contigs_2.fasta.gz
```

For four columns (with pre-annotated data):

```csv title="samplesheet.csv"
sample,fasta,protein,gbk
sample_1,/<path>/<to>/wastewater_metagenome_contigs_1.fasta.gz,/<path>/<to>/wastewater_metagenome_contigs_1.faa,/<path>/<to>/wastewater_metagenome_contigs_1.fasta.gbk
sample_2,/<path>/<to>/wastewater_metagenome_contigs_2.fasta.gz,/<path>/<to>/wastewater_metagenome_contigs_2.faa,/<path>/<to>/wastewater_metagenome_contigs_2.fasta.gbk
```

| Column    | Description                                                                                                                                                                                                           |
| --------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This will be used to name all output files from the pipeline. Spaces in sample names are automatically converted to underscores (`_`).                                                            |
| `fasta`   | Path or URL to a gzipped or uncompressed FASTA file. Accepted file suffixes are: `.fasta`, `.fna`, or `.fa`, or any of these with `.gz`, e.g. `.fa.gz`.                                                               |
| `protein` | Optional path to a pre-generated amino acid FASTA file (`.faa`) containing protein annotations of `fasta`, optionally gzipped. Required to be supplied if `gbk` also given.                                           |
| `gbk`     | Optional path to a pre-generated annotation file in Genbank format (`.gbk`, or `.gbff`) format containing annotations information of `fasta`, optionally gzipped. Required to be supplied if `protein` is also given. |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

:::danger
We highly recommend performing quality control on input contigs before running the pipeline. You may not receive results for some tools if none of the contigs in a FASTA file reach certain thresholds. Check parameter documentation for relevant minimum contig parameters.

For example, ideally BGC screening requires contigs of at least 3,000 bp else downstream tools may crash.
:::

## Notes on screening tools and taxonomic classification

The implementation of some tools in the pipeline may have some particular behaviours that you should be aware of before you run the pipeline.

### MMseqs2

MMseqs2 is currently the only taxonomic classification tool used in the pipeline to assign a taxonomic lineage to the input contigs. The database used to assign the taxonomic lineage can either be:

- a custom based database created by the user using `mmseqs createdb` externally and beforehand. If this flag is assigned, this database takes precedence over the default database in ` mmseqs_databases_id`.

  ```bash
  --taxa_classification_mmseqs_databases_localpath 'path/to/mmsesqs_custom_database/dir'
  ```

- an MMseqs2 ready database. These databases were compiled by the developers of MMseqs2 and can be called using their labels. All available options can be found [here](https://github.com/soedinglab/MMseqs2/wiki#downloading-databases). Only use those databases that have taxonomy files available (i.e., Taxonomy == Yes). By default mmseqs2 in the pipeline uses '[Kalamari](https://github.com/lskatz/Kalamari)' and runs an aminoacid based alignment.

  ```bash
  --taxa_classification_mmseqs_databases_id 'Kalamari'
  ```

### antiSMASH

antiSMASH has a minimum contig parameter, in which only contigs of a certain length (or longer) will be screened. In cases where no hits are found in these, the tool ends successfully without hits. However if no contigs in an input file reach that minimum threshold, the tool will end with a 'failure' code, and cause the pipeline to crash.

When the annotation is run with Prokka, the resulting `.gbk` file passed to antiSMASH may produce the error `translation longer than location allows` and end the pipeline run. This Prokka bug has been reported before (see [discussion on GitHub](https://github.com/antismash/antismash/discussions/450)) and is not likely to be fixed soon.

:::warning
If antiSMASH is run for BGC detection, we recommend to **not** run Prokka for annotation but instead use the default annotation tool (Pyrodigal) or switch > to Prodigal, or (for bacteria only!) Bakta.
:::

## Databases and reference files

Various tools of nf-core/funcscan use databases and reference files to operate.

nf-core/funcscan offers the functionality to auto-download databases for you, and as these databases can be very large, and we suggest to store these files in a central place from where you can reuse them across pipeline runs.

If your infrastructure has internet access (particularly on compute nodes), we **highly recommend** allowing the pipeline to download these databases for you on a first run, saving these to your results directory with `--save_databases`, then moving these to a different location (in case you wish to delete the results directory of this first run). An exception to this is HMM files where no auto-downloading functionality is possible.

:::warning
We generally do not recommend downloading the databases yourself, as this can often be non-trivial to do!
:::

As a reference, we will describe below where and how you can obtain databases and reference files used for tools included in the pipeline.

### Bakta

nf-core/funcscan offers multiple tools for annotating input sequences. Bakta is a new tool touted as a bacteria-only successor to the well-established Prokka.

To supply the preferred Bakta database (and not have the pipeline download it for every new run), use the flag `--annotation_bakta_db_localpath`. The full or light Bakta database must be downloaded from the Bakta Zenodo archive, the link of which can be found on the [Bakta GitHub repository](https://github.com/oschwengers/bakta#database-download).

Once downloaded this must be untarred:

```bash
tar xvzf db.tar.gz
```

And then passed to the pipeline with:

```bash
--annotation_bakta_db_localpath /<path>/<to>/db/
```

:::info
The flag `--save_databases` saves the pipeline-downloaded databases in your results directory. You can then move these to a central cache directory of your choice for re-use in the future.
:::

### hmmsearch

nf-core/funcscan allows screening of sequences for functional genes associated with various natural product types via Hidden Markov Models (HMMs) using hmmsearch.

This requires supplying a list of HMM files ending in `.hmm`, that have models for the particular molecule(s) or BGCs you are interested in. You can download these files from places such as [PFAM](https://www.ebi.ac.uk/interpro/download/Pfam/) for antimicrobial peptides (AMP), or the antiSMASH GitHub repository for [biosynthetic gene cluster](https://github.com/antismash/antismash/tree/master/antismash/detection/hmm_detection/data) related HMMs, or create them yourself.

You should place all HMMs in a directory and supply them e.g. to AMP models:

```bash
--amp_hmmsearch_models '/<path>/<to>/<amp>/*.hmm'
```

### AMRFinderPlus

AMRFinderPlus relies on NCBI's curated Reference Gene Database and curated collection of Hidden Markov Models.

nf-core/funcscan will download this database for you, unless the path to a local version is given with:

```bash
--arg_amrfinderplus_db '/<path>/<to>/<amrfinderplus_db>/'
```

To obtain a local version of the database:

1. Install AMRFinderPlus from [bioconda](https://bioconda.github.io/recipes/ncbi-amrfinderplus/README.html?highlight=amrfinderplus). To ensure database compatibility, please use the same version as is used in your nf-core/funcscan release (check version in file `<installation>/<path>/funcscan/modules/nf-core/amrfinderplus/run/environment.yml`).
2. Run `amrfinder --update`, which will download the latest version of the AMRFinderPlus database to the default location (location of the AMRFinderPlus binaries/data). It creates a directory in the format YYYY-MM-DD.version (e.g., `<installation>/<path>/data/2024-01-31.1/`).

<details markdown="1">
<summary>AMR related files in the database folder</summary>

```console
<YYYY-MM-DD.v>/
├── AMR_CDS.*
├── AMR_DNA-Campylobacter.*
├── AMR_DNA-Clostridioides_difficile.*
├── AMR_DNA-Enterococcus_faecalis.*
├── AMR_DNA-Enterococcus_faecium.*
├── AMR_DNA-Escherichia.*
├── AMR_DNA-Neisseria.*
├── AMR_DNA-Salmonella.*
├── AMR_DNA-Staphylococcus_aureus.*
├── AMR_DNA-Streptococcus_pneumoniae.*
├── AMR.LIB.*
├── AMRProt.*
├── changes.txt
├── database_format_version.txt
├── fam.tab
├── taxgroup.tab
└── version.txt
```

</details>

:::info
The flag `--save_databases` saves the pipeline-downloaded databases in your results directory. You can then move these to a central cache directory of your choice for re-use in the future.
:::

### DeepARG

DeepARG requires a database of potential antimicrobial resistance gene sequences based on a consensus from UNIPROT, CARD, and ARDB.

nf-core/funcscan can download this database for you, however it is very slow and pipeline runtime will be improved if you download this separately and supply it to the pipeline.

You can either:

1. Install DeepARG from [bioconda](https://bioconda.github.io/recipes/deeparg/README.html?highlight=deeparg)
2. Run `deeparg download_data -o /<path>/<to>/<database_location>/`

Or download the files directly from

1. the [DeepARG FTP site](https://bench.cs.vt.edu/ftp/data/gustavo1/deeparg/database/)
2. the [DeepARG database Zenodo archive](https://zenodo.org/record/8280582)

Note that more recent database versions maybe available from the [ARGMiner service](https://bench.cs.vt.edu/argminer/#/home).

You can then supply the path to resulting database directory with:

```bash
--arg_deeparg_data '/<path>/<to>/<deeparg>/<db>/'
```

Note that if you supply your own database that is not downloaded by the pipeline, make sure to also supply `--arg_deeparg_data_version` along
with the version number so hAMRonization will correctly display the database version in the summary report.
:::info
The flag `--save_databases` saves the pipeline-downloaded databases in your results directory. You can then move these to a central cache directory of your choice for re-use in the future.
:::

### RGI

RGI requires the database CARD which can be downloaded by nf-core/funcscan or supplied by the user manually. To download and supply the database yourself, do:

1. Download [CARD](https://card.mcmaster.ca/latest/data)
2. Extract the archive.

You can then supply the path to resulting database directory with:

```bash
--arg_rgi_database '/<path>/<to>/<card>/'
```

:::info
The flag `--save_databases` saves the pipeline-downloaded databases in your results directory. You can then move these to a central cache directory of your choice for re-use in the future.
:::

### antiSMASH

antiSMASH requires several databases for the detection of potential biosynthetic gene cluster (BGC) sequences (ClusterBlast, MIBiG, Pfam, Resfams, TIGRFAMs).

nf-core/funcscan can download these databases for you, however this is very slow and pipeline runtime will be improved if you download them separately and supply them to the pipeline.

The same applies for the antiSMASH installation directory, which is also a required parameter for the pipeline when using containers, due to some slight incompatibility when using such engines.

To supply the database directories to the pipeline:

1. Install antiSMASH from [bioconda](https://bioconda.github.io/recipes/antismash-lite/README.html)
2. Run `download-antismash-databases`
3. You can then supply the paths to the resulting databases and the whole installation directory with:

```bash
--bgc_antismash_databases '/<path>/<to>/<antismash>/<db>/'
--bgc_antismash_installationdirectory '/<path>/<to>/<antismash>/<dir>/'
```

Note that the names of the supplied folders must differ from each other (e.g. `antismash_db` and `antismash_dir`). If they are not provided, the databases will be auto-downloaded upon each BGC screening run of the pipeline.

:::info
The flag `--save_databases` saves the pipeline-downloaded databases in your results directory. You can then move these to a central cache directory of your choice for re-use in the future.
:::

:::info
If installing with conda, the installation directory will be `lib/python3.10/site-packages/antismash` from the base directory of your conda install or conda environment directory.
:::

### DeepBGC

DeepBGC relies on trained models and Pfams to run its analysis. nf-core/funcscan will download these databases for you. If the flag `--save_databases` is set, the downloaded files will be stored in the output directory under `databases/deepbgc/`.

Alternatively, if you already downloaded the database locally with `deepbgc download`, you can indicate the path to the database folder with `--bgc_deepbgc_database <path>/<to>/<deepbgc_db>/`. The folder has to contain the subfolders as in the database folder downloaded by `deepbgc download`:

```console
deepbgc_db/
├── common
  └── Pfam-hmm-models*.hmm.*
└── <version-num>[0.1.0]
  ├── classifier
  | └── myClassifiers*.pkl
  └── detector
    └── myDetectors*.pkl
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/funcscan -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

## Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/funcscan
```

## Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/funcscan releases page](https://github.com/nf-core/funcscan/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
