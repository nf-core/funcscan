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

When switched on, all tools of the given workflow will be run by default. If you don't need specific tools, you can explicitly skip them. The exception is HMMsearch, which needs to be explicitly switched on and provided with HMM screening files (AMP and BGC workflows, see [parameter documentation](/funcscan/parameters)). For the taxonomic classification, MMseqs2 is currently the only tool implemented in the pipline.

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

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/funcscan -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

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

## Notes on screening tools, taxonomic and functional classifications

The implementation of some tools in the pipeline may have some particular behaviours that you should be aware of before you run the pipeline.

### MMseqs2

MMseqs2 is currently the only taxonomic classification tool used in the pipeline to assign a taxonomic lineage to the input contigs. The database used to assign the taxonomic lineage can either be:

- A custom based database created by the user using `mmseqs createdb` externally and beforehand. If this flag is assigned, this database takes precedence over the default database in `--mmseqs_db_id`.

  ```bash
  --taxa_classification_mmseqs_db '<path>/<to>/<mmsesqs_custom_database>/<directory>'
  ```

  The contents of the directory should have files such as `<dbname>.version` and `<dbname>.taxonomy` in the top level.

- An MMseqs2 ready database. These databases were compiled by the developers of MMseqs2 and can be called using their labels. All available options can be found [here](https://github.com/soedinglab/MMseqs2/wiki#downloading-databases). Only use those databases that have taxonomy files available (i.e., Taxonomy == Yes). By default mmseqs2 in the pipeline uses '[Kalamari](https://github.com/lskatz/Kalamari)', and runs an aminoacid based alignment. However, if the user requires a more comprehensive taxonomic classification, we recommend the use of [GTDB](https://gtdb.ecogenomic.org/), but for that please remember to increase the memory, CPU threads and time required for the process `MMSEQS_TAXONOMY`.

  ```bash
  --taxa_classification_mmseqs_db_id 'Kalamari'
  ```

### InterProScan

[InterProScan](https://github.com/ebi-pf-team/interproscan) is currently the only protein annotation tool that gives a snapshot of the protein families and domains for each coding region.
By giving `--run_protein_annotation_interproscan`, the [InterPro database](http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.67-99.0/) v5.67-99.0 is by default downloaded and prepared and the input sequences will be screened against the database.
You can skip database downloading by the pipeline on each run by manually downloading and extracting the files from any [InterPro version](http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/) and giving the resulting directory path to `--protein_annotation_interproscan_db`.

```bash
--function_interproscan_db 'path/to/InterPro_directory/'
```

:::info
By default the databases used to assign the nearest protein domain is set as `PANTHER,ProSiteProfiles,ProSitePatterns,Pfam`. An addition of other application to the list, does not guarantee that the results will be integrated correctly within `AMPcombi`.
:::

### antiSMASH

antiSMASH has a minimum contig parameter, in which only contigs of a certain length (or longer) will be screened. In cases where no hits are found in these, the tool ends successfully without hits. However if no contigs in an input file reach that minimum threshold, the tool will end with a 'failure' code, and cause the pipeline to crash.

When the annotation is run with Prokka, the resulting `.gbk` file passed to antiSMASH may produce the error `translation longer than location allows` and end the pipeline run. This Prokka bug has been reported before (see [discussion on GitHub](https://github.com/antismash/antismash/discussions/450)) and is not likely to be fixed soon.

:::warning
If antiSMASH is run for BGC detection, we recommend to **not** run Prokka for annotation but instead use the default annotation tool (Pyrodigal) or switch to Prodigal, or (for bacteria only!) Bakta.
:::

## Databases and reference files

Various tools of nf-core/funcscan use databases and reference files to operate.

nf-core/funcscan offers the functionality to auto-download databases for you, and as these databases can be very large, and we suggest to store these files in a central place from where you can reuse them across pipeline runs.

If your infrastructure has internet access (particularly on compute nodes), we **highly recommend** allowing the pipeline to download these databases for you on a first run, saving these to your results directory with `--save_db`, then moving these to a different location (in case you wish to delete the results directory of this first run). An exception to this is HMM files where no auto-downloading functionality is possible.

:::warning
We generally do not recommend downloading the databases yourself, as this can often be non-trivial to do!
:::

As a reference, we will describe below where and how you can obtain databases and reference files used for tools included in the pipeline.

### Bakta

nf-core/funcscan offers multiple tools for annotating input sequences. Bakta is a new tool touted as a bacteria-only successor to the well-established Prokka.

To supply the preferred Bakta database (and not have the pipeline download it for every new run), use the flag `--annotation_bakta_db`.
The full or light Bakta database must be downloaded from the Bakta Zenodo archive.

You can do this by installing via conda and using the dedicated command

```bash
conda create -n bakta -c bioconda bakta
conda activate bakta

bakta_db download --output <LOCATION_TO_STORE> --type <full|light>
```

Alternatively, you can manually download the files via the links which can be found on the [Bakta GitHub repository](https://github.com/oschwengers/bakta#database-download).

Once downloaded this must be untarred:

```bash
tar xvzf db.tar.gz
```

And then passed to the pipeline with:

```bash
--annotation_bakta_db /<path>/<to>/<db>/
```

The contents of the directory should have files such as `*.dmnd` in the top level.

:::info
The flag `--save_db` saves the pipeline-downloaded databases in your results directory. You can then move these to a central cache directory of your choice for re-use in the future.
:::

### hmmsearch

nf-core/funcscan allows screening of sequences for functional genes associated with various natural product types via Hidden Markov Models (HMMs) using hmmsearch.

This requires supplying a list of HMM files ending in `.hmm`, that have models for the particular molecule(s) or BGCs you are interested in.
You can download these files from places such as [PFAM](https://www.ebi.ac.uk/interpro/download/Pfam/) for antimicrobial peptides (AMP), or the antiSMASH GitHub repository for [biosynthetic gene cluster](https://github.com/antismash/antismash/tree/master/antismash/detection/hmm_detection/data) related HMMs, or create them yourself.

You should place all HMMs in a directory, supply them to the AMP or BGC workflow and switch hmmsearch on:

```bash
--amp_run_hmmsearch --amp_hmmsearch_models "/<path>/<to>/<amp>/*.hmm"
```

:::warning
Ensure to wrap this path in double quotes if using an asterisk, to ensure Nextflow (not your shell) parses the wildcard.
:::

### AMPcombi

For AMPcombi, nf-core/funcscan will by default download the most recent version of the [DRAMP](http://dramp.cpu-bioinfor.org/) database as a reference database, and modifies the files for aligning the AMP hits in the AMP workflow.

nf-core/funcscan currently provides a python3 helper script to do these steps.

```bash
mkdir -p ampcombi/amp_ref_database
cd ampcombi/
wget https://github.com/nf-core/funcscan/raw/<PIPELINE_VERSION>/bin/ampcombi_download.py
python3 ampcombi_download.py
```

In addition to [DRAMP](http://dramp.cpu-bioinfor.org/), two more reference databases can be used to classify the recovered AMPs in the AMP workflow; [APD](https://aps.unmc.edu/) and [UniRef100](https://academic.oup.com/bioinformatics/article/23/10/1282/197795). Only one database can be used at a time using `--amp_ampcombi_db database_name`.

However, the user can also supply their own custom AMP database by following the guidelines in [AMPcombi](https://ampcombi.readthedocs.io/en/main/).
This can then be passed to the pipeline with:

```bash
--amp_ampcombi_db_dir_path '/<path>/<to>/<ampcombi_database>
```

The contents of the directory should have files such as `*.fasta` and `*.tsv` in the top level; a fasta file and the corresponding table with structural, functional and (if reported) taxonomic classifications. AMPcombi will then generate the corresponding `mmseqs2` directory, in which all binary files are prepared for downstream alignment of the recovered AMPs with [MMseqs2](https://github.com/soedinglab/MMseqs2). These can also be provided by the user by setting up an mmseqs2 compatible database using `mmseqs createdb *.fasta` in a directory called `mmseqs2`. An example file structure for [DRAMP](http://dramp.cpu-bioinfor.org/) used as the reference database:

```bash
amp_DRAMP_database/
├── general_amps_2024_11_13.fasta
├── general_amps_2024_11_13.txt
└── mmseqs2
    ├── ref_DB
    ├── ref_DB.dbtype
    ├── ref_DB_h
    ├── ref_DB_h.dbtype
    ├── ref_DB_h.index
    ├── ref_DB.index
    ├── ref_DB.lookup
    └── ref_DB.source
```

:::note{.fa-whale}
For both [DRAMP](http://dramp.cpu-bioinfor.org/) and [APD](https://aps.unmc.edu/), AMPcombi removes entries that contains any non amino acid residues by default.
:::

:::warning
The pipeline will automatically run Pyrodigal instead of Prodigal if the parameters `--run_annotation_tool prodigal --run_amp_screening` are both provided.
This is due to an incompatibility issue of Prodigal's output `.gbk` file with multiple downstream tools.
:::

:::tip

- If `--run_protein_annotation_interproscan` is given, protein and domain classifications of the coding regions are generated and the output is then integrated into the `AMPcombi parsetables` resulting table for every sample and the complete summary files e.g., `Ampcombi_summary.tsv`.

- In some cases when the AMP and the taxonomic classification subworkflows are
turned on, it can happen that only summary files per sample are created in the
output folder with **NO** `Ampcombi_summary.tsv` and `Ampcombi_summary_cluster.
tsv` files with no taxonomic classifications merged. This can occur if some AMP
parameters are 'too strict' or only one AMP tool is run, which can lead to no AMP
hits found in any of the samples or in only one sample. Look out for `[nf-core/
funcscan] AMPCOMBI2: 0/1 file passed. Skipping AMPCOMBI2_COMPLETE,
AMPCOMBI2_CLUSTER, and TAXONOMY MERGING steps.`in the stdout or `.nextflow.log`
file. In that case we recommend to lower the amp threshold and run more than one
AMP prediction tool.
:::

### Abricate

The default ABRicate installation comes with a series of 'default' databases:

- NCBI AMRFinderPlus (`ncbi`)
- CARD (`card`)
- ResFinder (`resfinder`)
- ARG-ANNOT (`argannot`)
- MEGARES (`megares`)
- EcOH (`echo`)
- PlasmidFinder (`plasmidfinder`)
- VFDB (`vfdb`)
- Ecoli_VF (`ecoli_vf`)

Each can be specified by using the nf-core/funcscan flag, for example for card: `--arg_abricate_db_id card`.

ABRicate also allows you to download additional and/or use custom databases.
For both of these, you will need to have your own local installation of ABRicate.
You then can download/add the custom database to the local installation's database directory, and supply this directory to the pipeline with the flag `--arg_abricate_db`, in combination with the name of the new database to `--arg_abricate_db_id <db_name>`.

For example, if you want to use the `bacmet2` database that does not come with the default installation, you could do:

```bash
## Create conda environment
conda create -n abricate -c bioconda abricate
conda activate abricate

## Download the bacmet2 database
abricate-get_db --db bacmet2 ## the logging will tell you where the database is downloaded to, e.g. /home/<user>/bin/miniconda3/envs/abricate/db/bacmet2/sequences
```

The resulting directory and database name can be passed to the pipeline as follows

```bash
--arg_abricate_db /<path>/<to>/<abricate>/db/ --arg_abricate_db_id bacmet2
```

The contents of the directory should have a directory named with the database name in the top level (e.g. `bacmet2/`).

### AMRFinderPlus

AMRFinderPlus relies on NCBI's curated Reference Gene Database and curated collection of Hidden Markov Models.

nf-core/funcscan will download this database for you, unless the path to a local version is given with:

```bash
--arg_amrfinderplus_db '/<path>/<to>/<amrfinderplus_db>/latest'
```

You must give the `latest` directory to the pipeline, and the contents of the directory should include files such as `*.nbd`, `*.nhr`, `versions.txt` etc. in the top level.

To obtain a local version of the database:

1. Install AMRFinderPlus from [bioconda](https://bioconda.github.io/recipes/ncbi-amrfinderplus/README.html?highlight=amrfinderplus).
   To ensure database compatibility, please use the same version as is used in your nf-core/funcscan release (check version in file `<installation>/<path>/funcscan/modules/nf-core/amrfinderplus/run/environment.yml`).

```bash
conda create -n amrfinderplus -c bioconda ncbi-amrfinderplus=3.12.8
conda activate amrfinderplus
```

2. Run `amrfinder --update`, which will download the latest version of the AMRFinderPlus database to the default location (location of the AMRFinderPlus binaries/data).
   It creates a directory in the format YYYY-MM-DD.version (e.g., `<installation>/<path>/data/2024-01-31.1/`).

<details markdown="1">
<summary>AMR related files in the database folder</summary>

```tree
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
The flag `--save_db` saves the pipeline-downloaded databases in your results directory. You can then move these to a central cache directory of your choice for re-use in the future.
:::

### DeepARG

DeepARG requires a database of potential antimicrobial resistance gene sequences based on a consensus from UNIPROT, CARD, and ARDB.

nf-core/funcscan can download this database for you, however it is very slow and pipeline runtime will be improved if you download this separately and supply it to the pipeline.

You can either:

1. Install DeepARG from [bioconda](https://bioconda.github.io/recipes/deeparg/README.html?highlight=deeparg)

```bash
conda create -n deeparg -c bioconda deeparg
conda activate deeparg
```

2. Run `deeparg download_data -o /<path>/<to>/<database_location>/`

Or download the files directly from

1. the [DeepARG FTP site](https://bench.cs.vt.edu/ftp/data/gustavo1/deeparg/database/)
2. the [DeepARG database Zenodo archive](https://zenodo.org/record/8280582)

Note that more recent database versions maybe available from the [ARGMiner service](https://bench.cs.vt.edu/argminer/#/home).

You can then supply the path to resulting database directory with:

```bash
--arg_deeparg_db '/<path>/<to>/<deeparg>/<db>/'
```

The contents of the directory should include directories such as `database`, `model`, and files such as `deeparg.gz` etc. in the top level.

Note that if you supply your own database that is not downloaded by the pipeline, make sure to also supply `--arg_deeparg_db_version` along
with the version number so hAMRonization will correctly display the database version in the summary report.

:::info
The flag `--save_db` saves the pipeline-downloaded databases in your results directory.
You can then move these to a central cache directory of your choice for re-use in the future.
:::

### MMSeqs2

To download MMSeqs2 databases for taxonomic classification, you can install `mmseqs` via conda:

```bash
conda create -n mmseqs2 -c bioconda mmseqs2
conda activate mmseqs2
```

Then to download the database of your choice

```bash
mmseqs databases <DATABASE_NAME> <LOCATION_TO_STORE> tmp/
```

:::info
You may want to specify a different location for `tmp/`, we just borrowed here from the official `mmseqs` [documentation](https://github.com/soedinglab/mmseqs2/wiki#downloading-databases).
:::

### RGI

RGI requires the database CARD which can be downloaded by nf-core/funcscan or supplied by the user manually.
To download and supply the database yourself, do:

1. Download [CARD](https://card.mcmaster.ca/latest/data)

```bash
wget https://card.mcmaster.ca/latest/data
```

2. Extract the (`.tar.bz2`) archive.

```bash
tar -xjvf data
```

You can then supply the path to resulting database directory with:

```bash
--arg_rgi_db '/<path>/<to>/<card>/'
```

The contents of the directory should include files such as `card.json`, `aro_index.tsv`, `snps.txt` etc. in the top level.

:::info
The flag `--save_db` saves the pipeline-downloaded databases in your results directory.
You can then move these to a central cache directory of your choice for re-use in the future.
:::

### antiSMASH

antiSMASH requires several databases for the detection of potential biosynthetic gene cluster (BGC) sequences (ClusterBlast, MIBiG, Pfam, Resfams, TIGRFAMs).

nf-core/funcscan can download these databases for you, however this is very slow and pipeline runtime will be improved if you download them separately and supply them to the pipeline.

The same applies for the antiSMASH installation directory, which is also a required parameter for the pipeline when using containers, due to some slight incompatibility when using such engines.

To supply the database directories to the pipeline:

1. Install antiSMASH from [bioconda](https://bioconda.github.io/recipes/antismash-lite/README.html). To ensure database compatibility, please use the same version as is used in your nf-core/funcscan release (check version in file `<pipeline_installation>/<path>/funcscan/modules/nf-core/antismash/antismashlite/environment.yml`).

```bash
conda create -n antismash-lite -c bioconda antismash-lite
conda activate antismash-lite
```

2. Run the command `download-antismash-databases`. Use `--database-dir` to specify a new location.
3. You can then supply the paths to the resulting databases and the whole installation directory with:

```bash
--bgc_antismash_db '/<path>/<to>/<antismash>/<db>/'
--bgc_antismash_installdir '/<path>/<to>/<antismash>/<dir>/antismash'
```

Note that the names of the supplied folders must differ from each other (e.g. `antismash_db` and `antismash_dir`).
The contents of the database directory should include directories such as `as-js/`, `clusterblast/`, `clustercompare/` etc. in the top level.
The contents of the installation directory should include directories such as `common/` `config/` and files such as `custom_typing.py` `custom_typing.pyi` etc. in the top level.

:::info
If installing with conda, the installation directory will be `lib/python3.10/site-packages/antismash` from the base directory of your conda install or conda environment directory.
:::

Note that the names of the two required folders must differ from each other (i.e., the `--bgc_antismash_db` directory must not be called `antismash`).
If they are not provided, the databases will be auto-downloaded upon each BGC screening run of the pipeline.

:::info
The flag `--save_db` saves the pipeline-downloaded databases in your results directory. You can then move these to a central cache directory of your choice for re-use in the future.
:::

### DeepBGC

DeepBGC relies on trained models and Pfams to run its analysis.
nf-core/funcscan will download these databases for you. If the flag `--save_db` is set, the downloaded files will be stored in the output directory under `databases/deepbgc/`.

Alternatively, you can download the database locally with:

```bash
conda create -n deepbgc -c bioconda deepbgc
conda activate deepbgc
export DEEPBGC_DOWNLOADS_DIR=<PREFERRED_CACHE_DIRECTORY>
deepbgc download
```

You can then indicate the path to the database folder in the pipeline with `--bgc_deepbgc_db <path>/<to>/<deepbgc_db>/`.
The contents of the database directory should include directories such as `common`, `0.1.0` in the top level:

```console
deepbgc_db/
├── common
A diifferent version of the database can be supplied to the pipeline by passing the InterProScan database directory to `--protein_annotation_interproscan_db path/to/interproscan_db/`. The directory can be created following with:
└── <version-num>[0.1.0]
  ├── classifier
  | └── myClassifiers*.pkl
  └── detector
    └── myDetectors*.pkl
```

### InterProScan

[InterProScan](https://github.com/ebi-pf-team/interproscan) is used to provide more information about the proteins annotated on the contigs. By default, turning on this subworkflow with `--run_protein_annotation_interproscan` will download and unzip the (as of now) latest [InterPro database](http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.67-99.0/) v5.67-99.0. The database downloaded can be saved in the output directory `<output_directors>/databases/interproscan/*` if the `--save_db` is turned on. Note: the download can take upto 4 hours depending on teh bandwidth.

A different version of the database can be supplied to the pipeline by passing the InterProScan database directory to `--protein_annotation_interproscan_db path/to/downloaded-untarred-interproscan_db-dir/`. The directory can be created following:

```
curl -L https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.67-99.0/interproscan-5.67-99.0-64-bit.tar.gz -o interproscan_db/interproscan-5.67-99.0-64-bit.tar.gz
tar -xzf interproscan_db/interproscan-5.67-99.0-64-bit.tar.gz -C interproscan_db/

```

## Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/funcscan
```

## Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/funcscan releases page](https://github.com/nf-core/funcscan/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

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

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

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
