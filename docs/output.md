# nf-core/funcscan: Output

## Introduction

The output of nf-core/funcscan provides the output directories from each tool applied, as well as a summary of tool outputs for each of the functional groups: antibiotic resistance genes ([AMRfinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/),[DeepARG](https://bitbucket.org/gusphdproj/deeparg-ss/src/master/), [fargene](https://github.com/fannyhb/fargene), [RGI](https://card.mcmaster.ca/analyze/rgi)), antimicrobial peptides ([macrel](https://github.com/BigDataBiology/macrel), [amplify](https://github.com/bcgsc/AMPlify), [ampir](https://ampir.marine-omics.net/), [hmmsearch](http://hmmer.org)), biosynthetic gene clusters ([antiSMASH](https://docs.antismash.secondarymetabolites.org)) and functional annotation ([prokka](https://github.com/tseemann/prokka)) and ([prodigal](https://github.com/hyattpd/Prodigal)).

Furthermore, for reproducibility, versions of all software used in the run is presented in a [MultiQC](http://multiqc.info) report.

The directories listed below will be created in the results directory (specified by the `--outdir` flag) after the pipeline has finished. All paths are relative to this top-level output directory. The default directory sturcture of nf-core/funcscan is:

<!--
```bash
outdir/
#├── acep/
#├── ai4amp/
├── antismash/
├── amplify/
├── ampir/
#├── combiamp/
├── deeparg/
#├── ensembleamppred/
├── fargene/
├── hamronizer/
├── hmmsearch/
├── macrel/
├── multiqc/
#├── neubi/
├── pipeline_info/
├── prodigal/
├── prokka/
#└── rgi/
work/
```
-->

```bash
outdir/
├── antismash/
├── amplify/
├── ampir/
├── amrfinder/
├── deeparg/
├── fargene/
├── hamronizer/
├── hmmsearch/
├── macrel/
├── multiqc/
├── pipeline_info/
├── prodigal/
├── prokka/
└── rgi/
work/
```

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

Antimicrobial Resistance Genes (ARGs):

- [AMRfinderPlus](#amrfinderplus) - antimicrobial resistance gene detection, using NCBI’s curated Reference Gene Database and curated collection of Hidden Markov Models
- [DeepARG](#deeparg) - antimicrobial resistance gene detection, using a deep learning model
- [fARGene](#fargene) - antimicrobial resistance gene detection, using a deep learning model
- [rgi](#rgi) - antimicrobial resistance gene detection, based on alignment to the CARD database

Antimicrobial Peptides (AMPs) and peptide annotation:

- [Prodigal](#prodigal) - for open reading frame annotation
- [Prokka](#prokka) - (optional: alternative to prodigal) open reading frame and functional protein annotation
  <!--* [acep](#acep) - antimicrobial peptide detection-->
  <!--* [ai4amp](#ai4amp) - antimicrobial peptide detection-->
- [ampir](#ampir) - antimicrobial peptide detection
- [amplify](#amplify) - antimicrobial peptide detection
  <!--* [EnsembleAMPPred](#ensembleamppred) - antimicrobial peptide detection-->
- [Macrel](#macrel) - antimicrobial peptide detection
<!--* [neubi](#neubi) - antimicrobial peptide detection-->

Biosynthetic Gene Clusters (BGCs):

- [hmmsearch](#hmmsearch) - biosynthetic gene cluster detection
- [antiSMASH](#antismash) - biosythedic gene cluster detection

Output Summaries:

- [MultiQC](#multiqc) - Report of all software and versions used in the pipeline.
- [hAMRonization](#hamronization) - summarizes resistance gene output from various detection tools
  <!--* [combiAMP](#combiamp) - summarizes antimicrobial peptide detection output-->
  <!--* [comBGC](#combgc) - PRELIMINARY TOOL NAME - summarizes biosynthetic gene cluster detection output-->
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing raw parsed data used for MultiQC report rendering.
  - `multiqc_plots/`: directory containing any static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is used in nf-core/funcscan to report the versions of all software used in the given pipeline run. This allows for reproducible analysis and transparency in method reporting in publications.

### hAMRonization

<details markdown="1">
<summary>Output files</summary>

- `hamronization/` one of the following:
  - `hamronization_combined_report.json`: summarized output in .json format
  - `hamronization_combined_report.tsv`: summarized output in .tsv format
  - `hamronization_combined_report.html`: interactive output in .html format

</details>

[hAMRonization](https://github.com/pha4ge/hAMRonization) summarizes the outputs of the **antimicrobial resistance gene** detection tools (DeepARG, fARGene) into a single unified format. It supports a variety of summary options including an interactive summary.

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

### AMRfinderPlus

<details markdown="1">
<summary>Output files</summary>

- `amrfinderplus/`
  - `db/`: contains the database downloaded with `amrfinderplus update`
  - `*.tsv`: contains the search results

</details>

[AMRfinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) relies on NCBI’s curated Reference Gene Database and curated collection of Hidden Markov Models. It identifies antimicrobial resistance genes, resistance-associated point mutations, and select other classes of genes using protein annotations and/or assembled nucleotide sequence.

The `*.tsv` output files contain the following fields:

- Protein identifier
- Contig id
- Start
- Stop
- Strand
- Gene symbol
- Sequence name
- Scope
- Element type
- Element subtype
- Class
- Subclass
- Method
- Target length
- Reference sequence length
- % Coverage of reference sequence
- % Identity to reference sequence
- Alignment length
- Accession of closest sequence
- Name of closest sequence
- HMM id
- HMM description

### DeepARG

<details markdown="1">
<summary>Output files</summary>

- `deeparg/`
  - `db/`: contains diamond, data, database and model information
  - `predict/`:
    - `*.align.daa*`: Diamond alignment output.
    - `*.align.daa.tsv`: Diamond alignment output as .tsv.
    - `*.mapping.ARG`: contains the sequences with a probability >= --prob (0.8 default).
    - `*.mapping.potential.ARG`: contains the sequences with a probability < --prob (0.8 default).

</details>

[deepARG](https://bitbucket.org/gusphdproj/deeparg-ss/src/master/) uses deep learning to characterize and annotate antibiotic resistance genes in metagenomes. It is composed of two models for two types of input: short sequence reads and gene-like sequences. In this pipeline we use the `ls` model, which is suitable for annotating full sequence genes and to discover novel antibiotic resistance genes from assembled samples. The tool `Diamond` is used as an aligner.

The `*.ARG` output files contain the following fields:

- ARG_NAME
- QUERY_START
- QUERY_END
- QUERY_ID
- PREDICTED_ARG_CLASS
- BEST_HIT_FROM_DATABASE
- PREDICTION_PROBABILITY
- ALIGNMENT_BESTHIT_IDENTITY (%)
- ALIGNMENT_BESTHIT_LENGTH
- ALIGNMENT_BESTHIT_BITSCORE
- ALIGNMENT_BESTHIT_EVALUE
- COUNTS

### fARGene

<details markdown="1">
<summary>Output files</summary>

- `fargene/`
  - `fargene_analysis.log`: Contains the output that Fargene produced during its run
  - `<sample_name>/`:
    - `hmmsearchresults/`: Contains the output from hmmsearch.
    - `predictedGenes/`:
      - `*-filtered.fasta`: nucleotide sequences of predicted ARGs.
      - `*-filtered-peptides.fasta`: aminoacid sequences of predicted ARGs.
    - `results_summary.txt`: Text summary of run results, listing predicted genes and ORFs for each input file.
    - `tmpdir/`: Contains temporary output files and fasta files.

</details>

[fARGene](https://github.com/fannyhb/fargene) (Fragmented Antibiotic Resistance Gene Identifier) is a tool that takes either fragmented metagenomic data or longer sequences as input and predicts and delivers full-length antiobiotic resistance genes as output. The tool includes developed and optimized models for a number or resistance gene types. The model to be used has to be specified for each run (`--hmm-model`). Models available are:

- `class_a`: Class A beta-lactamases
- `class_b_1_2`: Subclass B1 and B2 beta-lactamases
- `class_b3`: Subclass B3 beta-lactamases
- `class_d_1`: Class C beta-lactamases
- `class_d_2`: Class D beta-lactamases
- `qnr`: Quinolone resistance genes
- `tet_efflux`, `tet_rpg`, `tet_enzyme`: Tetracycline resistance genes

### RGI

<details markdown="1">
<summary>Output files</summary>

- `rgi/`
  - `<samplename>.json`: hit results in json format
  - `<samplename>.txt`: hit results table separated by '#'
  - `<samplename>.fasta.temp.*.json`: four temporary json files where '\*' stands for 'homolog', 'overexpression', 'predictedGenes' and 'predictedGenes.protein'.

</details>

[RGI](https://github.com/arpcard/rgi) (Resistance Gene Identifier) predicts resistome(s) from protein or nucleotide data based on homology and SNP models. It uses reference data from the Comprehensive Antibiotic Resistance Database (CARD).

### Prodigal

<details markdown="1">
<summary>Output files</summary>

- `prokka/`
  - `<samplename>/`:
    - `*.gff`: Annotation in GFF3 format, containing both sequences and annotations
    - `*.fna`: Nucleotide FASTA file of the input contig sequences.
    - `*.faa`: Protein FASTA file of the translated CDS sequences.
    - `*_all.txt`: Text file containing all_gene_annotations.

</details>

[Prodigal](https://github.com/hyattpd/Prodigal) is an alternative for prokka that does whole genome annotation to identify CDS in a set of genomic DNA sequences. It can be applied to annotate bacterial, archaeal and viral genomes.

### Prokka

<details markdown="1">
<summary>Output files</summary>

- `prokka/`
  - `<samplename>/`:
    - `*.gff`: annotation in GFF3 format, containing both sequences and annotations
    - `*.gbk`: standard Genbank file derived from the master .gff.
    - `*.fna`: Nucleotide FASTA file of the input contig sequences.
    - `*.faa`: Protein FASTA file of the translated CDS sequences.
    - `*.ffn`: Nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA).
    - `*.sqn`: An ASN1 format "Sequin" file for submission to Genbank.
    - `*.fsa`: Nucleotide FASTA file of the input contig sequences, used by "tbl2asn" to create the .sqn file.
    - `*.tbl`: Feature Table file, used by "tbl2asn" to create the .sqn file.
    - `*.err`: Unacceptable annotations - the NCBI discrepancy report.
    - `*.log`: Contains all the output that Prokka produced during its run.
    - `*.txt`: Statistics relating to the annotated features found.
    - `*.tsv`: ab-separated file of all features.

</details>

[Prokka](https://github.com/tseemann/prokka) does whole genome annotation to identify features of interest in a set of genomic DNA sequences, and labelling them with useful information. It can be applied to annotate bacterial, archaeal and viral genomes.

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

### Ampir

<details markdown="1">
<summary>Output files</summary>

- `ampir/`
  - `<samplename>.ampir.faa`: predicted AMP sequences in FAA format
  - `<samplename>.ampir.tsv`: predicted AMP metadata in TSV format, contains contig name, sequence and probability score

</details>

[ampir](https://github.com/Legana/ampir) (antimicrobial peptide prediction in r) package was designed to predict antimicrobial peptides (AMPs) from any given size protein dataset. ampir uses a supervised statistical machine learning approach to predict AMPs. It incorporates two support vector machine classification models, “precursor” and “mature” that have been trained on publicly available antimicrobial peptide data.

### AMPlify

<details markdown="1">
<summary>Output files</summary>

-_`amplify/` -_ `*_results.tsv`: contig amino-acid sequences with prediction result (AMP or non-AMP) and information on sequence length, charge, probability score, AMPlify log-scaled score)

</details>

[AMPlify](https://github.com/bcgsc/AMPlify) is an attentive deep learning model for antimicrobial peptide prediction. It takes in annotated contigs (.faa) and classifies them as either AMP or non-AMP.

<!--### Ensemble-AMPPred

<details markdown="1">
<summary>Output files</summary>

* `ensembleamppred/`
    * `output1`: xxx
    * `output2/`: xxx

</details>

[Ensemble-AMPPred](no link to source code found ...) xxx tool description here xxx
-->

### Macrel

<details markdown="1">
<summary>Output files</summary>

- `macrel_contigs/`
  - `*.smorfs.faa.gz`: A zipped fasta file containing aminoacid sequences of small peptides (<100 aa, small open reading frames) showing the general gene prediction information in the contigs.
  - `*.all_orfs.faa.gz`: A zipped fasta file containing amino acid sequences showing the general gene prediction information in the contigs.
  - `prediction.gz`: A zipped file, with all predicted amps in a table format.
  - `*.md`: A readme file containing tool specific information (e.g. citations, details about the output, etc.).
  - `*_log.txt`: A log file containing the information pertaining to the run.

</details>

[Macrel](https://github.com/BigDataBiology/macrel) is a tool that mines antimicrobial peptides (AMPs) from (meta)genomes by predicting peptides from genomes (provided as contigs) and outputs all the predicted anti-microbial peptides found.

<!--### NeuBI

<details markdown="1">
<summary>Output files</summary>

* `neubi/`
    * `output1`: xxx
    * `output2/`: xxx

</details>

[NeuBI](https://github.com/nafizh/NeuBI) (Neural Bacteriocin Identifier) is a recurrent neural network based software to predict bacteriocins from protein sequences. Unlike traditional alignment based approaches such as BLAST or HMMER used by BAGEL or BACTIBASE, this is an alignment free approach towards finding novel bacteriocins.
-->

### hmmsearch

<details markdown="1">
<summary>Output files</summary>

- `hmmersearch/`
  - `*.txt.gz`: Human readable output summarizing hmmsearch results.
  - `*.sto.gz`: Optional multiple sequence alignment (MSA) in Stockholm format.
  - `*.tbl.gz`: Optional tabular (space-delimited) summary of per-target output.
  - `*.domtbl.gz`: Optional tabular (space-delimited) summary of per-domain output.

</details>

[HMMER/hmmsearch](http://hmmer.org) is used for searching sequence databases for sequence homologs, and for making sequence alignments. It implements methods using probabilistic models called profile hidden Markov models (profile HMMs). `hmmsearch` is used to search one or more profiles against a sequence database.

### antiSMASH

<details markdown="1">
<summary>Output files</summary>

- `antismash/` most important output files:
  - `knownclusterblast/`
    - `*_c*.txt`: Tables with MIBiG hits
  - `clusterblastoutput.txt`: Raw BLAST output of known clusters previously predicted by antiSMASH using the built-in ClusterBlast algorithm
  - `knownclusterblastoutput.txt`: Raw BLAST output of known clusters of the MIBiG database.
  - `*region*.gbk`: Nucleotide sequence + annotations in GenBank file format; one file per antiSMASH hit.

</details>

[antiSMASH](https://docs.antismash.secondarymetabolites.org) (antibiotics & Secondary Metabolite Analysis Shell) is a tool for rapid genome-wide identification, annotation and analysis of secondary metabolite biosynthesis gene clusters in bacterial and fungal genomes. It identifies biosynthetic loci covering the whole range of known secondary metabolite compound classes and aligns the identified regions at the gene cluster level to their nearest relatives from a database containing all other known gene clusters. It integrates or cross-links all previously available secondary-metabolite specific gene analysis methods in one interactive view.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
