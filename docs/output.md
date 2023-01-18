# nf-core/funcscan: Output

## Introduction

The output of nf-core/funcscan provides reports for each of the functional groups:

- antibiotic resistance genes ([ABRicate](https://github.com/tseemann/abricate), [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder), [DeepARG](https://bitbucket.org/gusphdproj/deeparg-ss/src/master), [fARGene](https://github.com/fannyhb/fargene), [RGI](https://card.mcmaster.ca/analyze/rgi))
- antimicrobial peptides ([macrel](https://github.com/BigDataBiology/macrel), [amplify](https://github.com/bcgsc/AMPlify), [ampir](https://ampir.marine-omics.net), [hmmsearch](http://hmmer.org))
- biosynthetic gene clusters ([antiSMASH](https://docs.antismash.secondarymetabolites.org), [deepBGC](https://github.com/Merck/deepbgc), [GECCO](https://gecco.embl.de), [hmmsearch](http://hmmer.org))

Additional to summary reports, the output directories from all applied tools are provided. This includes the functional annotation output from either [prokka](https://github.com/tseemann/prokka), [prodigal](https://github.com/hyattpd/Prodigal), or [Bakta](https://github.com/oschwengers/bakta) if the `--save_annotations` flag was set.

Similarly, all downloaded databases are saved (i.e. from [antiSMASH](https://docs.antismash.secondarymetabolites.org), [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder), [Bakta](https://github.com/oschwengers/bakta), and/or [DeepARG](https://bitbucket.org/gusphdproj/deeparg-ss/src/master)) into the output directory `<outdir>/downloads/` if the `--save_databases` flag was set.

Furthermore, for reproducibility, versions of all software used in the run is presented in a [MultiQC](http://multiqc.info) report.

The directories listed below will be created in the results directory (specified by the `--outdir` flag) after the pipeline has finished. All paths are relative to this top-level output directory. The default directory structure of nf-core/funcscan is:

<!--
```console
results/
├── annotation/
|   ├── prodigal/
|   └── prokka/
├── amp/
#|   ├── acep/
#|   ├── ai4amp/
|   ├── ampir/
|   ├── amplify/
#|   ├── ensembleamppred/
|   ├── hmmsearch/
|   ├── macrel/
#|   └── neubi/ (https://github.com/nafizh/NeuBI)
├── arg/
|   ├── abricate/
|   ├── amrfinderplus/
|   ├── deeparg/
|   ├── fargene/
|   ├── hamronization/
|   └── rgi/
├── bgc/
|   ├── antismash/
|   ├── gecco/
|   └── hmmsearch/
├── reports/
|   ├── AMPcombi/
#|   ├── BGCcombi/
|   └── hamronization_summarize/
├── databases/
├── multiqc/
└── pipeline_info/
work/
```
-->

```console
results/
├── annotation/
|   ├── prodigal/
|   ├── prokka/
|   └── bakta/
├── amp/
|   ├── ampir/
|   ├── amplify/
|   ├── hmmsearch/
|   └── macrel/
├── arg/
|   ├── abricate/
|   ├── amrfinderplus/
|   ├── deeparg/
|   ├── fargene/
|   ├── hamronization/
|   └── rgi/
├── bgc/
|   ├── antismash/
|   ├── deepbgc/
|   ├── gecco/
|   └── hmmsearch/
├── reports/
|   ├── ampcombi/
|   ├── comBGC/
|   └── hamronization_summarize/
├── databases/
├── multiqc/
└── pipeline_info/
work/
```

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data through the following steps:

DNA sequence annotation

- [Prodigal](#prodigal) – for open reading frame annotation.
- [Prokka](#prokka) – (optional: alternative to prodigal) open reading frame and functional protein annotation.
- [Bakta](#bakta) – (optional: alternative to prokka) open reading frame and functional protein annotation.

Antimicrobial Resistance Genes (ARGs):

- [ABRicate](#abricate) – antimicrobial resistance gene detection, based on alignment to one of several databases.
- [AMRFinderPlus](#amrfinderplus) – antimicrobial resistance gene detection, using NCBI’s curated Reference Gene Database and curated collection of Hidden Markov Models.
- [DeepARG](#deeparg) – antimicrobial resistance gene detection, using a deep learning model.
- [fARGene](#fargene) – antimicrobial resistance gene detection, using a deep learning model.
- [RGI](#rgi) – antimicrobial resistance gene detection, based on alignment to the CARD database.

Antimicrobial Peptides (AMPs):

  <!--* [acep](#acep) – antimicrobial peptide detection-->
  <!--* [ai4amp](#ai4amp) – antimicrobial peptide detection-->

- [ampir](#ampir) – antimicrobial peptide detection, based on a supervised statistical machine learning approach.
- [amplify](#amplify) – antimicrobial peptide detection, using a deep learning model.
  <!--* [EnsembleAMPPred](#ensembleamppred) – antimicrobial peptide detection-->
- [hmmsearch](#hmmsearch) – antimicrobial peptide detection, based on hidden Markov models.
- [Macrel](#macrel) – antimicrobial peptide detection, using a machine learning approach.
<!--* [neubi](#neubi) – antimicrobial peptide detection-->

Biosynthetic Gene Clusters (BGCs):

- [antiSMASH](#antismash) – biosythetic gene cluster detection.
- [deepBGC](#deepbgc) - biosynthetic gene cluster detection, using a deep learning model.
- [GECCO](#gecco) – biosynthetic gene cluster detection, using using Conditional Random Fields (CRFs).
- [hmmsearch](#hmmsearch) – biosynthetic gene cluster detection, based on hidden Markov models.

Output Summaries:

- [MultiQC](#multiqc) – report of all software and versions used in the pipeline.
- [hAMRonization](#hamronization) – summary of resistance gene output from various detection tools.
- [AMPcombi](#ampcombi) – summary of antimicrobial peptide output from various detection tools.
- [comBGC](#combgc) – summary of biosynthetic gene cluster detection output
- [Pipeline information](#pipeline-information) – report metrics generated during the workflow execution.

## Tool details

### Prodigal

<details markdown="1">
<summary>Output files</summary>

- `prodigal/`
  - `<samplename>/`:
    - `*.gff`: annotation in GFF3 format, containing both sequences and annotations
    - `*.fna`: nucleotide FASTA file of the input contig sequences
    - `*.faa`: protein FASTA file of the translated CDS sequences
    - `*.gbk`: annotation in GBK format, containing both sequences and annotations

</details>

[Prodigal](https://github.com/hyattpd/Prodigal) is an alternative for prokka that does whole genome annotation to identify CDS in a set of genomic DNA sequences. It can be applied to annotate bacterial, archaeal and viral genomes.

### Prokka

<details markdown="1">
<summary>Output files</summary>

- `prokka/`
  - `<samplename>/`
    - `*.gff`: annotation in GFF3 format, containing both sequences and annotations
    - `*.gbk`: standard Genbank file derived from the master .gff
    - `*.fna`: nucleotide FASTA file of the input contig sequences
    - `*.faa`: protein FASTA file of the translated CDS sequences
    - `*.ffn`: nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA)
    - `*.sqn`: an ASN1 format "Sequin" file for submission to Genbank
    - `*.fsa`: nucleotide FASTA file of the input contig sequences, used by "tbl2asn" to create the .sqn file
    - `*.tbl`: feature Table file, used by "tbl2asn" to create the .sqn file
    - `*.err`: unacceptable annotations - the NCBI discrepancy report
    - `*.log`: logging output that Prokka produced during its run
    - `*.txt`: statistics relating to the annotated features found
    - `*.tsv`: tab-separated file of all features

</details>

[Prokka](https://github.com/tseemann/prokka) does whole genome annotation to identify features of interest in a set of genomic DNA sequences, and labelling them with useful information. It can be applied to annotate bacterial, archaeal and viral genomes.

### Bakta

<details markdown="1">
<summary>Output files</summary>

- `bakta/`
  - `<samplename>.gff3`: annotations & sequences in GFF3 format
  - `<samplename>.gbff`: annotations & sequences in (multi) GenBank format
  - `<samplename>.ffn`: feature nucleotide sequences as FASTA
  - `<samplename>.fna`: replicon/contig DNA sequences as FASTA
  - `<samplename>.embl`: annotations & sequences in (multi) EMBL format
  - `<samplename>.faa`: CDS/sORF amino acid sequences as FASTA
  - `<samplename>_hypothetical.faa`: further information on hypothetical protein CDS as simple human readble tab separated values
  - `<samplename>_hypothetical.tsv`: hypothetical protein CDS amino acid sequences as FASTA
  - `<samplename>.tsv`: annotations as simple human readble TSV
  - `<samplename>.txt`: summary in TXT format

> Descriptions directly from the [Bakta documentation](https://github.com/oschwengers/bakta#output).

</details>

[Bakta](https://github.com/oschwengers/bakta) is a tool for the rapid & standardized annotation of bacterial genomes and plasmids from both isolates and MAGs. It provides dbxref-rich, sORF-including and taxon-independent annotations in machine-readable JSON & bioinformatics standard file formats for automated downstream analysis.

### ampir

<details markdown="1">
<summary>Output files</summary>

- `ampir/`
  - `<samplename>.ampir.faa`: predicted AMP sequences in FAA format
  - `<samplename>.ampir.tsv`: predicted AMP metadata in TSV format, contains contig name, sequence and probability score

</details>

[ampir](https://github.com/Legana/ampir) (**a**nti**m**icrobial **p**eptide **p**rediction **i**n **r**) was designed to predict antimicrobial peptides (AMPs) from any given size protein dataset. ampir uses a supervised statistical machine learning approach to predict AMPs. It incorporates two support vector machine classification models, “precursor” and “mature” that have been trained on publicly available antimicrobial peptide data.

### AMPlify

<details markdown="1">
<summary>Output files</summary>

- `amplify/`
  - `*_results.tsv`: contig amino-acid sequences with prediction result (AMP or non-AMP) and information on sequence length, charge, probability score, AMPlify log-scaled score)

</details>

[AMPlify](https://github.com/bcgsc/AMPlify) is an attentive deep learning model for antimicrobial peptide prediction. It takes in annotated contigs (as protein sequences) and classifies them as either AMP or non-AMP.

### hmmsearch

<details markdown="1">
<summary>Output files</summary>

- `hmmersearch/`
  - `*.txt.gz`: human readable output summarizing hmmsearch results
  - `*.sto.gz`: optional multiple sequence alignment (MSA) in Stockholm format
  - `*.tbl.gz`: optional tabular (space-delimited) summary of per-target output
  - `*.domtbl.gz`: optional tabular (space-delimited) summary of per-domain output

</details>

[HMMER/hmmsearch](http://hmmer.org) is used for searching sequence databases for sequence homologs, and for making sequence alignments. It implements methods using probabilistic models called profile hidden Markov models (profile HMMs). `hmmsearch` is used to search one or more profiles against a sequence database.

### Macrel

<details markdown="1">
<summary>Output files</summary>

- `macrel_contigs/`
  - `*.smorfs.faa.gz`: zipped fasta file containing amino acid sequences of small peptides (<100 aa, small open reading frames) showing the general gene prediction information in the contigs
  - `*.all_orfs.faa.gz`: zipped fasta file containing amino acid sequences showing the general gene prediction information in the contigs
  - `prediction.gz`: zipped file, with all predicted and non-predicted amps in a table format
  - `*.md`: readme file containing tool specific information (e.g. citations, details about the output, etc.)
  - `*_log.txt`: log file containing the information pertaining to the run

</details>

[Macrel](https://github.com/BigDataBiology/macrel) is a tool that mines antimicrobial peptides (AMPs) from (meta)genomes by predicting peptides from genomes (provided as contigs) and outputs all the predicted anti-microbial peptides found.

### AMPcombi

<details markdown="1">
<summary>Output files</summary>

- `ampcombi/`
  - `ampcombi_complete_summary.csv`: summarized output from all AMP workflow tools (except hmmer_hmmsearch) in csv
  - `ampcombi.log`: a log file generated by ampcombi
  - `*_ampcombi.csv`: summarized output in csv for each sample
  - `*_amp.faa*`: fasta file containing the amino acid sequences for all AMP hits for each sample
  - `*_diamond_matches.txt*`: alignment file generated fby DIAMOND for each sample

</details>

[AMPcombi](https://github.com/Darcy220606/AMPcombi) summarizes the results of **antimicrobial peptide (AMP)** prediction tools (AMPIR, AMPLIFY, MACREL, and other non nf-core tools) into a single table and aligns the hits against a reference AMP database for functional and taxonomic classification.

### ABRicate

<details markdown="1">
<summary>Output files</summary>

- `abricate/`
  - `*.{csv,tsv}`: search results in tabular format

</details>

[ABRicate](https://github.com/tseemann/abricate) screens contigs for antimicrobial resistance or virulence genes. It comes bundled with multiple databases: NCBI, CARD, ARG-ANNOT, Resfinder, MEGARES, EcOH, PlasmidFinder, Ecoli_VF and VFDB.

### AMRFinderPlus

<details markdown="1">
<summary>Output files</summary>

- `amrfinderplus/`
  - `*.tsv`: search results in tabular format

</details>

[AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder) relies on NCBI’s curated Reference Gene Database and curated collection of Hidden Markov Models. It identifies antimicrobial resistance genes, resistance-associated point mutations, and select other classes of genes using protein annotations and/or assembled nucleotide sequence.

### DeepARG

<details markdown="1">
<summary>Output files</summary>

- `deeparg/`
  - `*.align.daa*`: Diamond alignment output
  - `*.align.daa.tsv`: Diamond alignment output as .tsv
  - `*.mapping.ARG`: ARG predictions with a probability >= --prob (0.8 default).
  - `*.mapping.potential.ARG`: ARG predictions with a probability < --prob (0.8 default)

</details>

[deepARG](https://bitbucket.org/gusphdproj/deeparg-ss/src/master) uses deep learning to characterize and annotate antibiotic resistance genes in metagenomes. It is composed of two models for two types of input: short sequence reads and gene-like sequences. In this pipeline we use the `ls` model, which is suitable for annotating full sequence genes and to discover novel antibiotic resistance genes from assembled samples. The tool `Diamond` is used as an aligner.

### fARGene

<details markdown="1">
<summary>Output files</summary>

- `fargene/`
  - `fargene_analysis.log`: logging output that Fargene produced during its run
  - `<sample_name>/`:
    - `hmmsearchresults/`: output from hmmsearch (only if `--arg_fargene_savetmpfiles` supplied)
    - `predictedGenes/`:
      - `*-filtered.fasta`: nucleotide sequences of predicted ARGs
      - `*-filtered-peptides.fasta`: aminoacid sequences of predicted ARGs
    - `results_summary.txt`: text summary of results, listing predicted genes and ORFs for each input file
    - `tmpdir/`: temporary output files and fasta files (only if `--arg_fargene_savetmpfiles` supplied)

</details>

[fARGene](https://github.com/fannyhb/fargene) (**F**ragmented **A**ntibiotic **R**esistance **G**ene Identifier) is a tool that takes either fragmented metagenomic data or longer sequences as input and predicts and delivers full-length antibiotic resistance genes as output. The tool includes developed and optimized models for a number of resistance gene types. By default the pipeline runs all models. If only a sub-list or single model is to be used this can be specified with the `--hmm-model` flag. Available models are:

- `class_a`: class A beta-lactamases
- `class_b_1_2`: subclass B1 and B2 beta-lactamases
- `class_b3`: subclass B3 beta-lactamases
- `class_d_1`: class C beta-lactamases
- `class_d_2`: class D beta-lactamases
- `qnr`: quinolone resistance genes
- `tet_efflux`, `tet_rpg`, `tet_enzyme`: tetracycline resistance genes

### RGI

<details markdown="1">
<summary>Output files</summary>

- `rgi/`
  - `<samplename>.json`: hit results in json format
  - `<samplename>.txt`: hit results table separated by '#'
  - `temp/`:
    - `<samplename>.fasta.temp.*.json`: temporary json files, '\*' stands for 'homolog', 'overexpression', 'predictedGenes' and 'predictedGenes.protein' (only if `--arg_rgi_savetmpfiles` supplied).

</details>

[RGI](https://github.com/arpcard/rgi) (**R**esistance **G**ene **I**dentifier) predicts resistome(s) from protein or nucleotide data based on homology and SNP models. It uses reference data from the Comprehensive Antibiotic Resistance Database (CARD).

### hAMRonization

<details markdown="1">
<summary>Output files</summary>

- `hamronization/` one of the following:
  - `hamronization_combined_report.json`: summarized output in .json format
  - `hamronization_combined_report.tsv`: summarized output in .tsv format
  - `hamronization_combined_report.html`: interactive output in .html format

</details>

[hAMRonization](https://github.com/pha4ge/hAMRonization) summarizes the outputs of the **antimicrobial resistance gene** detection tools (ABRicate, AMRFinderPlus, DeepARG, fARGene, RGI) into a single unified format. It supports a variety of summary options including an interactive summary.

### antiSMASH

<details markdown="1">
<summary>Output files</summary>

- `antismash/`
  - `knownclusterblast/`
    - `*_c*.txt`: tables with MIBiG hits
  - `clusterblastoutput.txt`: raw BLAST output of known clusters previously predicted by antiSMASH using the built-in ClusterBlast algorithm
  - `knownclusterblastoutput.txt`: raw BLAST output of known clusters of the MIBiG database.
  - `*region*.gbk`: nucleotide sequence + annotations in GenBank file format; one file per antiSMASH hit.

</details>

[antiSMASH](https://docs.antismash.secondarymetabolites.org) (**anti**biotics & **S**econdary **M**etabolite **A**nalysis **Sh**ell) is a tool for rapid genome-wide identification, annotation and analysis of secondary metabolite biosynthesis gene clusters in bacterial and fungal genomes. It identifies biosynthetic loci covering the whole range of known secondary metabolite compound classes and aligns the identified regions at the gene cluster level to their nearest relatives from a database containing all other known gene clusters. It integrates or cross-links all previously available secondary-metabolite specific gene analysis methods in one interactive view.

### deepBGC

<details markdown="1">
<summary>Output files</summary>

- `deepbgc/`
  - `README.txt`: Summary of output files generated
  - `LOG.txt`: Log output of DeepBGC
  - `*.antismash.json`: AntiSMASH JSON file for sideloading
  - `*.bgc.gbk`: Sequences and features of all detected BGCs in GenBank format
  - `*.bgc.tsv`: Table of detected BGCs and their properties
  - `*.full.gbk`: Fully annotated input sequence with proteins, Pfam domains (PFAM_domain features) and BGCs (cluster features)
  - `*.pfam.tsv`: Table of Pfam domains (pfam_id) from given sequence (sequence_id) in genomic order, with BGC detection scores
  - `evaluation/`
    - `*.bgc.png`: Detected BGCs plotted by their nucleotide coordinates
    - `*.pr.png`: Precision-Recall curve based on predicted per-Pfam BGC scores
    - `*.roc.png`: ROC curve based on predicted per-Pfam BGC scores
    - `*.score.png`: BGC detection scores of each Pfam domain in genomic order

</details>

[deepBGC](https://github.com/Merck/deepbgc) detects BGCs in bacterial and fungal genomes using deep learning. DeepBGC employs a Bidirectional Long Short-Term Memory Recurrent Neural Network and a word2vec-like vector embedding of Pfam protein domains. Product class and activity of detected BGCs is predicted using a Random Forest classifier.

### GECCO

<details markdown="1">
<summary>Output files</summary>

- `gecco/`
  - `*.genes.tsv/`: TSV file containing detected/predicted genes with BGC probability scores
  - `*.features.tsv`: TSV file containing identified domains
  - `*.clusters.tsv`: TSV file containing coordinates of predicted clusters and BGC types
  - `*_cluster_*.gbk`: GenBank file (if clusters were found) containing sequence with annotations; one file per GECCO hit
  - `*.json`: antiSMASH v6 sideload JSON file (if `--antismash-sideload`) supplied

</details>

[GECCO](https://gecco.embl.de) is a fast and scalable method for identifying putative novel Biosynthetic Gene Clusters (BGCs) in genomic and metagenomic data using Conditional Random Fields (CRFs).

### comBGC
<details markdown="1">
<summary>Output files</summary>

- `comBGC/`
  - `bgc_summary.csv`: summary of the output of the pipelines' BGC detection tools used.

</details>

**comBGC** is a local module which summarizes the results of the **Biosynthetic Gene Cluster (BGC)** prediction tools (antiSMASH, deepBGC, GECCO) used in the pipeline into one comprehensive summary with standardized headers.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser
  - `multiqc_data/`: directory containing raw parsed data used for MultiQC report rendering
  - `multiqc_plots/`: directory containing any static images from the report in various formats

</details>

[MultiQC](http://multiqc.info) is used in nf-core/funcscan to report the versions of all software used in the given pipeline run. This allows for reproducible analysis and transparency in method reporting in publications.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

<!--### Ensemble-AMPPred

<details markdown="1">
<summary>Output files</summary>

* `ensembleamppred/`
    * `output1`: xxx
    * `output2/`: xxx

</details>

[Ensemble-AMPPred](no link to source code found ...) xxx tool description here xxx
-->

<!--### CombiAMP

<details markdown="1">
<summary>Output files</summary>

* `combiamp/`
    * `output1`: xxx
    * `output2/`: xxx

</details>

[CombiAMP](https://link-to-tool-page.org) xxx tool description here xxx SUMMARY of AMP tools' output
-->

<!--### ComBGC

<details markdown="1">
<summary>Output files</summary>

* `combiamp/`
    * `output1`: xxx
    * `output2/`: xxx

</details>

[ComBGC](https://link-to-tool-page.org) xxx tool description here xxx SUMMARY of BGC tools' output
-->

<!--### Acep

<details markdown="1">
<summary>Output files</summary>

* `acep/`
    * `output1`: xxx
    * `output2/`: xxx

</details>

[Acep](no page with source code found ...) xxx tool description here xxx
-->

<!--### AI4AMP

<details markdown="1">
<summary>Output files</summary>

* `ai4amp/`
    * `output1`: xxx
    * `output2/`: xxx

</details>

[AI4AMP](https://github.com/LinTzuTang/AI4AMP_predictor) is a sequence-based antimicrobial peptides (AMP) predictor based on PC6 protein encoding method and deep learning.
-->

<!--### NeuBI

<details markdown="1">
<summary>Output files</summary>

* `neubi/`
    * `output1`: xxx
    * `output2/`: xxx

</details>

[NeuBI](https://github.com/nafizh/NeuBI) (Neural Bacteriocin Identifier) is a recurrent neural network based software to predict bacteriocins from protein sequences. Unlike traditional alignment based approaches such as BLAST or HMMER used by BAGEL or BACTIBASE, this is an alignment free approach towards finding novel bacteriocins.
-->
