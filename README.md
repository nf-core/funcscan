# ![nf-core/funscan](docs/images/nf-core-funcscan_logo_flat_light.png#gh-light-mode-only) ![nf-core/funscan](docs/images/nf-core-funcscan_logo_flat_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/nf-core/funcscan/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/funcscan/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/funcscan/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/funcscan/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?logo=Amazon%20AWS)](https://nf-co.re/funcscan/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/funcscan)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23funcscan-4A154B?logo=slack)](https://nfcore.slack.com/channels/funcscan)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/funcscan** is a bioinformatics best-practice analysis pipeline for the screening of functional components of nucleotide sequences such as assembled contigs. This includes mining for antimicrobial peptides, antibiotic resistance genes and biosynthetic gene clusters.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/funcscan/results).

## Pipeline summary

1. Annotates prokaryotic input assembled contigs ([`Prodigal`](https://github.com/hyattpd/Prodigal) or [`PROKKA`](https://github.com/tseemann/prokka))
2. Screens contigs for antimicrobial peptide-like sequences with: [ampir](https://cran.r-project.org/web/packages/ampir/index.html), [Macrel](https://github.com/BigDataBiology/macrel), [HMMER](http://hmmer.org/), [AMPlify](https://github.com/bcgsc/AMPlify)
3. Screens contigs for antibiotic resistant gene-like sequences with:. [fARGene](https://github.com/fannyhb/fargene), [RGI](https://card.mcmaster.ca/analyze/rgi), [DEEPARG](https://bench.cs.vt.edu/deeparg)
4. Screens contigs for biosynthetic gene cluster-like sequences with: <!-- [antiSMASH](https://antismash.secondarymetabolites.org/#!/start), [gecco](https://gecco.embl.de/), [HMMER](http://hmmer.org/) -->
5. Software version reporting with ([`MultiQC`](http://multiqc.info/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow run nf-core/funcscan -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```console
   nextflow run nf-core/funcscan --input samplesheet.csv --outdir <OUTDIR> -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --run_amp_screening --run_arg_screening
   ```

## Documentation

The nf-core/funcscan pipeline comes with documentation about the pipeline [usage](https://nf-co.re/funcscan/usage), [parameters](https://nf-co.re/funcscan/parameters) and [output](https://nf-co.re/funcscan/output).

## Credits

nf-core/funcscan was originally written by Jasmin Frangenberg, Anan Ibrahim, Louisa Perelo, Moritz E. Beber, James A. Fellows Yates.

We thank the following people for their extensive assistance in the development of this pipeline:

Rosa Herbst, Martin Klapper.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#funcscan` channel](https://nfcore.slack.com/channels/funcscan) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/funcscan for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
