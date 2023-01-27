# nf-core/funcscan: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/funcscan/usage](https://nf-co.re/funcscan/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

nf-core/funcscan is a pipeline for efficient and parallelised screening of long nucleotide sequences such as contigs for possible natural products sequences. It targets antimicrobial peptides, antimicrobial resistance genes, and biosynthetic clusters.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/funcscan --input samplesheet.csv --outdir <OUTDIR> -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing temporary files required for the run
<OUTDIR>            # Final results (location specified with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs
```

To run any of the three screening workflows, switch them on:

- `--run_amp_screening`
- `--run_arg_screening`
- `--run_bgc_screening`

When switched on, all tools of the given workflow will be run by default. If you don't need specific tools, you can explicitly skip them.

Example: You want to run AMP and ARG screening but you don't need the DeepARG tool of the ARG workflow and the Macrel tool of the AMP workflow. Your command would be:

```bash
nextflow run nf-core/funcscan --input samplesheet.csv --outdir <OUTDIR> -profile docker --run_arg_screening --arg_skip_deeparg --run_amp_screening --arg_skip_macrel
```

## Samplesheet input

nf-core/funcscan takes FASTA files as input, typically contigs or whole genome sequences. To supply these to the pipeline, you will need to create a samplesheet with information about the samples you would like to analyses. Use this parameter to specify its location.

```bash
--input '[path to samplesheet file]'
```

The input samplesheet has to be a comma-separated file (`.csv`) with 2 columns (`sample`, and `fasta`), and a header row as shown in the examples below.

```bash
sample,fasta
sample_1,https://raw.githubusercontent.com/nf-core/test-datasets/funcscan/wastewater_metagenome_contigs_1.fasta.gz
sample_2,https://raw.githubusercontent.com/nf-core/test-datasets/funcscan/wastewater_metagenome_contigs_2.fasta.gz
```

| Column   | Description                                                                                                                                                |
| -------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample` | Custom sample name. This will be used to name all output files from the pipeline. Spaces in sample names are automatically converted to underscores (`_`). |
| `fasta`  | Full path to a gzipped or uncompressed FASTA file. Accepted file suffixes are: `.fasta`, `.fna`, or `.fa`, or any of these with `.gz`, e.g. `.fa.gz`.      |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

> ⚠️ We highly recommend performing quality control on input contigs before running the pipeline. You may not receive results for some tools if none of the contigs in a FASTA file reach certain thresholds. Check parameter documentation for relevent minimum contig parameters.

## Notes on screening tools

The implementation of some tools in the pipeline may have some particular behaviours that you should be aware of before you run the pipeline.

### antiSMASH

AntiSMASH has a minimum contig parameter, in which only contigs of a certain length (or longer) will be screened. In cases where no hits are found in these, then the tool ends successfully. However if no contigs in an input file reach that minimum threshold, the tool will end with a 'failure' code, and cause the pipeline to crash.

To prevent entire pipeline failures to due a single 'bad sample', nf-core/funcscan will filter out any input sample in which none of the contigs reach the minimum contig length specified with `--bgc_antismash_sampleminlength`.

> ⚠️ If a sample does not reach this contig length threshold, you will recieve a warning in your console and `.nextflow.log` file, and no result files will exist for this sample in your results directory.

## Databases and reference files

nf-core/funcscan utilises various tools that use databases and reference files to generate results. While nf-core/funcscan offers in some cases functionality to autodownload databases for you, these databases can be very large, and it is more efficient to store these files in a central place from where they can be reused across pipeline runs.

Here we will describe where you can obtain databases and reference files for tools included in the pipeline.

### Bakta

nf-core/funcscan offers multiple tools for annotating input sequences.

Bakta is a new tool touted as a bacteria-only successor to the well-established Prokka.

Bakta requires a database, supplied to the pipeline with `--annotation_bakta_db`.
This must be downloaded from the Bakta Zenodo archive, the link of which can be found on the [Bakta github repository](https://github.com/oschwengers/bakta).

Once downloaded this must be untarred

```bash
tar xvzf db.tar.gz
```

And then passed to the pipeline with

```bash
--annotation_bakta_db /<path>/<to>/db/
```

### hmmsearch

nf-core/funcscan allows screening of sequences for various natural product types via hidden markov models with HMMSearch.

This requires supplying a list of HMM model files ending in `.hmm`, that have models for the particular molecule(s) you are interested in.

You can download such files from places such as [https://pfam.xfam.org/](https://pfam.xfam.org/) for antimicrobial peptides (AMP), or the antiSMASH github repository for [biosynthetic gene cluster](https://github.com/antismash/antismash/tree/master/antismash/detection/hmm_detection/data) related HMM models

You should place all of these in a directory and supply them e.g. to AMP models

```bash
--amp_hmmsearch_models '/<path>/<to>/<amp>/*.hmm'
```

### ARMfinderPlus

AMRFinderPlus relies on NCBI’s curated Reference Gene Database and curated collection of Hidden Markov Models.

nf-core/funcscan will download this database for you, unless the path to a local version is given with:

```bash
--arg_amrfinderplus_db '/<path>/<to>/<amrfinderplus_db>/'
```

You can either:

1. Install AMRFinderPlus from [bioconda](https://bioconda.github.io/recipes/ncbi-amrfinderplus/README.html?highlight=amrfinderplus)
2. Run `amrfinder --update`, which will download the latest version of the AMRFinderPlus database to the default location (location of the AMRFinderPlus binaries/data). It creates a directory under data in the format YYYY-MM-DD.version (e.g., 2019-03-06.1).

Or

1. Download the files directly form the [NCBI FTP site](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/)

The downloaded database folder contains the AMR related files:

```console
<YYYY-MM-DD.v>[2022-08-09.1]/
├── AMR.LIB.*
├── AMRProt.*
├── AMR_CDS.*
├── AMR_DNA-Campylobacter.*
├── AMR_DNA-Clostridioides_difficilie.*
├── AMR_DNA-Enterococcus_faecalis.*
├── AMR_DNA-Enterococcus_faecium.*
├── AMR_DNA-Escherichia.*
├── AMR_DNA-Neisseria.*
├── AMR_DNA-Salmonella.*
├── AMR_DNA-Staphylococcus_aureus.*
├──AMR_DNA-Streptococcus_pneumoniae.*
├── changes.txt
├── database_format_version.txt
├── fam.tab
├── taxgroup.tab
└── version.txt
```

The path to the resulting database directory can be supplied to the pipeline as described above.

### DeepARG

DeepARG requires a database of potential antimicrobial resistence gene sequences, based on a consensus from UNIPROT, CARD and ARDB.

> ⚠️ As of January 2023 the DeepARG database server is down! We have deactivated automated database downloading in the pipeline - you can try to download your own copy using the instructions below!

nf-core/funcscan can download this database for you, however it is very slow and pipeline runtime will be improved if you download this separately and supply it to the pipeline.

You can either:

1. Install DeepARG from [bioconda](https://bioconda.github.io/recipes/deeparg/README.html?highlight=deeparg)
2. Run `deeparg download_data -o /<path>/<to>/<database_location>/`

Or

1. Download the files directly from the [DeepARG FTP site](https://bench.cs.vt.edu/ftp/data/gustavo1/deeparg/database/)

Note that more recent database versions maybe available from the [ARGMiner service](https://bench.cs.vt.edu/argminer/#/home).

You can then supply the path to resulting database directory with:

```bash
--arg_deeparg_data '/<path>/<to>/<deeparg>/<db>/'
```

Note that if you supply your own database that is not downloaded by the pipeline, make sure to also supply `--arg_deeparg_data_version` along
with the version number so hAMRonization will correctly display the database version in the summary report.

### AntiSMASH

AntiSMASH requires several databases of potential biosynthetic gene cluster (BGC) sequences (ClusterBlast, MIBiG, Pfam, Resfams, TIGRFAMs).

nf-core/funcscan can download these databases for you, however this is very slow and pipeline runtime will be improved if you download them separately and supply them to the pipeline.
The same applies for the antiSMASH installation directory, which is also a required parameter for the pipeline, due to some slight incompatibility when using containers.

To supply the database directories to the pipeline:

1. Install antiSMASH from [bioconda](https://bioconda.github.io/recipes/antismash-lite/README.html)
2. Run `download-antismash-databases`
3. You can then supply the paths to the resulting databases and the whole installation directory with:

```bash
--bgc_antismash_databases '/<path>/<to>/<antismash>/<db>/'
--bgc_antismash_installationdirectory '/<path>/<to>/<antismash>/<dir>/'
```

Note that the names of the supplied folders must differ from each other (e.g. `antismash_db` and `antismash_dir`). If they are not provided, the databases will be auto-downloaded upon each BGC screening run of the pipeline.

Hint: The flag `--save_databases` saves the pipeline-downloaded databases in your results directory. You can then move these to a central cache directory of your choice for re-use in the future.

> If installing with conda, the installation directory will be `lib/python3.8/site-packages/antismash` from the base directory of your conda install or conda environment directory.

### DeepBGC

DeepBGC relies on trained models and Pfam database to run its analysis. nf-core/funcscan will download these databases for you. If the flag `--save_databases` is set, the downloaded files will be stored in the output directory under `databases/deepbgc/`.

Alternatively, if you already downloaded the database locally with `deepbgc download`, you can indicate the path to the database folder with `--bgc_deepbgc_database path/to/deepbgc_db/`. The folder has to contain the subfolders as in the database folder downloaded by `deepbgc download`:

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

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/funcscan
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/funcscan releases page](https://github.com/nf-core/funcscan/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

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
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

#### For beginners

A first step to bypass this error, you could try to increase the amount of CPUs, memory, and time for the whole pipeline. Therefor you can try to increase the resource for the parameters `--max_cpus`, `--max_memory`, and `--max_time`. Based on the error above, you have to increase the amount of memory. Therefore you can go to the [parameter documentation of rnaseq](https://nf-co.re/rnaseq/3.9/parameters) and scroll down to the `show hidden parameter` button to get the default value for `--max_memory`. In this case 128GB, you than can try to run your pipeline again with `--max_memory 200GB -resume` to skip all process, that were already calculated. If you can not increase the resource of the complete pipeline, you can try to adapt the resource for a single process as mentioned below.

#### Advanced option on process level

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/star/align/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers (advanced users)

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

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
