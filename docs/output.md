# nf-core/funcscan: Output

## Introduction

The output of nf-core/funcscan provides reports for each of the functional groups:

- antibiotic resistance genes (tools: [ABRicate](https://github.com/tseemann/abricate), [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder), [DeepARG](https://bitbucket.org/gusphdproj/deeparg-ss/src/master), [fARGene](https://github.com/fannyhb/fargene), [RGI](https://card.mcmaster.ca/analyze/rgi) summarized by [hAMRonization](https://github.com/pha4ge/hAMRonization))
- antimicrobial peptides (tools: [macrel](https://github.com/BigDataBiology/macrel), [amplify](https://github.com/bcgsc/AMPlify), [ampir](https://ampir.marine-omics.net), [hmmsearch](http://hmmer.org) summarized by [AMPcombi](https://github.com/Darcy220606/AMPcombi))
- biosynthetic gene clusters (tools: [antiSMASH](https://docs.antismash.secondarymetabolites.org), [deepBGC](https://github.com/Merck/deepbgc), [GECCO](https://gecco.embl.de), [hmmsearch](http://hmmer.org) summarized by [comBGC](#combgc))

As a general workflow, we recommend to first look at the to summary reports ([ARGs](#hamronization), [AMPs](#ampcombi), [BGCs](#combgc)), to get a general overview of what hits have been found across all the tools of each functional group. After which, you can explore the specific output directories  of each tool to get more detailed information about each results. The tool-specific output directories, also includes the functional annotation output from either [prokka](https://github.com/tseemann/prokka), [prodigal](https://github.com/hyattpd/Prodigal), or [Bakta](https://github.com/oschwengers/bakta) if the `--save_annotations` flag was set.

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
- [comBGC](#combgc) – summary of biosynthetic gene cluster output from various detection tools.
- [Pipeline information](#pipeline-information) – report metrics generated during the workflow execution.

## Tool details

### Annotation tools:

[Prodigal](#prodigal), [Prokka](#prokka), [Bakta](#bakta)

#### Prodigal

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

#### Prokka

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

#### Bakta

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

### AMP detection tools:

[ampir](#ampir), [AMPlify](#amplify), [hmmsearch](#hmmsearch), [Macrel](#macrel), summary: [AMPcombi](#ampcombi)

#### ampir

<details markdown="1">
<summary>Output files</summary>

- `ampir/`
  - `<samplename>.ampir.faa`: predicted AMP sequences in FAA format
  - `<samplename>.ampir.tsv`: predicted AMP metadata in TSV format, contains contig name, sequence and probability score

</details>

[ampir](https://github.com/Legana/ampir) (**a**nti**m**icrobial **p**eptide **p**rediction **i**n **r**) was designed to predict antimicrobial peptides (AMPs) from any given size protein dataset. ampir uses a supervised statistical machine learning approach to predict AMPs. It incorporates two support vector machine classification models, “precursor” and “mature” that have been trained on publicly available antimicrobial peptide data.

#### AMPlify

<details markdown="1">
<summary>Output files</summary>

- `amplify/`
  - `*_results.tsv`: contig amino-acid sequences with prediction result (AMP or non-AMP) and information on sequence length, charge, probability score, AMPlify log-scaled score)

</details>

[AMPlify](https://github.com/bcgsc/AMPlify) is an attentive deep learning model for antimicrobial peptide prediction. It takes in annotated contigs (as protein sequences) and classifies them as either AMP or non-AMP.

#### hmmsearch

<details markdown="1">
<summary>Output files</summary>

- `hmmersearch/`
  - `*.txt.gz`: human readable output summarizing hmmsearch results
  - `*.sto.gz`: optional multiple sequence alignment (MSA) in Stockholm format
  - `*.tbl.gz`: optional tabular (space-delimited) summary of per-target output
  - `*.domtbl.gz`: optional tabular (space-delimited) summary of per-domain output

</details>

[HMMER/hmmsearch](http://hmmer.org) is used for searching sequence databases for sequence homologs, and for making sequence alignments. It implements methods using probabilistic models called profile hidden Markov models (profile HMMs). `hmmsearch` is used to search one or more profiles against a sequence database.

#### Macrel

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

#### AMPcombi

<details markdown="1">
<summary>Output files</summary>

- `ampcombi/`
  - `ampcombi_complete_summary.csv`: summarized output from all AMP workflow tools (except hmmer_hmmsearch) in csv
  - `ampcombi.log`: a log file generated by ampcombi
  - `*_ampcombi.csv`: summarized output in csv for each sample
  - `*_amp.faa*`: fasta file containing the amino acid sequences for all AMP hits for each sample
  - `*_diamond_matches.txt*`: alignment file generated fby DIAMOND for each sample

</details>

<details markdown="1">
<summary>AMP summary table headers</summary>

| Table column | Description |
|---|---|
| `name` | the name of the sample. |
| `contig_id` | the contig header. |
| `prob_macrel` | the probability associated with the AMP prediction using `MACREL`
| `prob_neubi` | the probability associated with the AMP prediction using `NEUBI`
| `prob_ampir` | the probability associated with the AMP prediction using `AMPIR`
| `prob_amplify` | the probability associated with the AMP prediction using `AMPLIFY`
| `evalue_hmmer` | the expected number of false positives (nonhomologous sequences) with a  similar of higher score. This stands for how significant te hit is, the lower the evalue, the
more significant the hit. |
| `aa_sequence` | the amino acid sequence that forms part of teh contig and is AMP encoding. |
| `target_id` | the [DRAMP](http://dramp.cpu-bioinfor.org/) ID within the database found to be similar to the predicted AMP by `DIAMOND` alignment. |
| `pident` | the percentage identity of amino acid residues that fully aligned between the `DRAMP` sequence and the predicted AMP sequence. |
| `evalue` | the number of alignments of similar or better qualities that can be expected when searching a database of similar size with a random sequence distribution. This is generated by `DIAMOND` alignments using the [DRAMP](http://dramp.cpu-bioinfor.org/) AMP database. The lower the value the more significant that the hit is positive. An e-value of < 0.001 means that the this hit will be found by chance once per 1,0000 queries. |
| `Sequence` | the sequence corresponding to the `DRAMP` ID found to be similar to the predicted AMP sequence. |
| `Sequence_length` | the number of amino acid residues in the `DRAMP` sequence. |
| `Name` | the full name of the peptide copied from teh database it was uploaded to. |
| `Swiss_Prot_Entry` | the entry name of the peptide within the [UniProtKB/Swiss-Prot](https://www.uniprot.org/help/entry_name) database. |
| `Family` | the name of teh family, group or class of AMPs this peptide belongs to. For e.g. bacteriocins. |
| `Gene` | the name of the gene (if available in the database) that encodes the peptide. |
| `Source` | the name of the source organism (if available in the database) from which the peptide was extracted from. |
| `Activity` | the activity of the peptide. For example, whether it is Antimicrobial, Antibacterial, Anti-Gram+, Anti-Gram-, Insecticidal or Antifungal. |
| `Protein_existence` | the status of the peptide. For example, whether it is only a homology, protein level, predicted or transcript level. |
| `Structure` | specify what type of structure the peptide has. For e.g. alpha helix, bridge, etc. |
| `Structure_Description` | specify any further description of the structure if available. |
| `PDB_ID` | the ID of an equivalent peptide found in the protein data bank [PDB](https://www.rcsb.org/docs/general-help/organization-of-3d-structures-in-the-protein-data-bank). |
| `Comments` | provide any further details found in the databse regarding the peptide. |
| `Target_Organism` | the name of the target organsim to which the peptide is effective against. |
| `Hemolytic_activity` | specify the type of hemolytic activity if any. |
| `Linear/Cyclic/Branched` | specify whether the hit is a linear, cyclic or branched peptide. |
| `N-terminal_Modification` | specify whether it contains N-terminal_Modification. |
| `C-termina_Modification` | specify whether it contains C-terminal_Modification. |
| `Other_Modifications` | specify whether there are any other modifications found in the peptide structure. |
| `Stereochemistry`— the type of peptide stereochemistry if available. |
| `Cytotoxicity` | the mechansim of cytotoxicity of the peptide if available. |
| `Binding_Target` | the binding target of the peptide. For e.g. lipid, cell membran or chitin binding. |
| `Pubmed_ID` | the Pubmed ID if a publication is associated with the peptide. |
| `Reference` | the citation of the associated publication if available. |
| `Author` | the author names associated with the publication or who have uploaded the peptide. |
| `Title` | the publication title if available. |

</details>

[AMPcombi](https://github.com/Darcy220606/AMPcombi) summarizes the results of **antimicrobial peptide (AMP)** prediction tools (AMPIR, AMPLIFY, MACREL, and other non nf-core tools) into a single table and aligns the hits against a reference AMP database for functional and taxonomic classification.

### ARG detection tools:

[ABRicate](#abricate), [AMRFinderPlus](#amrfinderplus), [DeepARG](#deeparg), [fARGene](#fargene), [RGI](#rgi), summary: [hAMRonization](#hamronization)

#### ABRicate

<details markdown="1">
<summary>Output files</summary>

- `abricate/`
  - `*.{csv,tsv}`: search results in tabular format

</details>

[ABRicate](https://github.com/tseemann/abricate) screens contigs for antimicrobial resistance or virulence genes. It comes bundled with multiple databases: NCBI, CARD, ARG-ANNOT, Resfinder, MEGARES, EcOH, PlasmidFinder, Ecoli_VF and VFDB.

#### AMRFinderPlus

<details markdown="1">
<summary>Output files</summary>

- `amrfinderplus/`
  - `*.tsv`: search results in tabular format

</details>

[AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder) relies on NCBI’s curated Reference Gene Database and curated collection of Hidden Markov Models. It identifies antimicrobial resistance genes, resistance-associated point mutations, and select other classes of genes using protein annotations and/or assembled nucleotide sequence.

#### DeepARG

<details markdown="1">
<summary>Output files</summary>

- `deeparg/`
  - `*.align.daa*`: Diamond alignment output
  - `*.align.daa.tsv`: Diamond alignment output as .tsv
  - `*.mapping.ARG`: ARG predictions with a probability >= --prob (0.8 default).
  - `*.mapping.potential.ARG`: ARG predictions with a probability < --prob (0.8 default)

</details>

[deepARG](https://bitbucket.org/gusphdproj/deeparg-ss/src/master) uses deep learning to characterize and annotate antibiotic resistance genes in metagenomes. It is composed of two models for two types of input: short sequence reads and gene-like sequences. In this pipeline we use the `ls` model, which is suitable for annotating full sequence genes and to discover novel antibiotic resistance genes from assembled samples. The tool `Diamond` is used as an aligner.

#### fARGene

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

#### RGI

<details markdown="1">
<summary>Output files</summary>

- `rgi/`
  - `<samplename>.json`: hit results in json format
  - `<samplename>.txt`: hit results table separated by '#'
  - `temp/`:
    - `<samplename>.fasta.temp.*.json`: temporary json files, '\*' stands for 'homolog', 'overexpression', 'predictedGenes' and 'predictedGenes.protein' (only if `--arg_rgi_savetmpfiles` supplied).

</details>

[RGI](https://github.com/arpcard/rgi) (**R**esistance **G**ene **I**dentifier) predicts resistome(s) from protein or nucleotide data based on homology and SNP models. It uses reference data from the Comprehensive Antibiotic Resistance Database (CARD).

#### hAMRonization

<details markdown="1">
<summary>Output files</summary>

- `hamronization/` one of the following:
  - `hamronization_combined_report.json`: summarized output in .json format
  - `hamronization_combined_report.tsv`: summarized output in .tsv format
  - `hamronization_combined_report.html`: interactive output in .html format

</details>
<details markdown="1">
<summary>ARG summary table headers</summary>

| Table column | Description |
 |---|---|
 | Required | |
 | `input_file_name` | The name of the file containing the sequence data to be analysed |
 | `gene_symbol` | The short name of a gene; a single word that does not contain white space characters. It is typically derived from the gene name |
 | `gene_name` | The name of a gene |
 | `reference_database_name` | An identifier of a biological or bioinformatics database |
 | `reference_database_version` | The version of the database containing the reference sequences used for analysis |
 | `reference_accession` | An identifier that specifies an individual sequence record in a public sequence repository |
 | `analysis_software_name` | A name of a computer package, application, method or function used for the analysis of data |
 | `analysis_software_version` | The version of software used to analyze data |
 | `genetic_variation_type` | The class of genetic variation detected |
 | Optional | |
 | `antimicrobial_agent` | A substance that kills or slows the growth of microorganisms, including bacteria, viruses, fungi and protozoans |
 | `coverage_percentage` | The percentage of the reference sequence covered by the sequence of interest |
 | `coverage_depth` | The average number of reads representing a given nucleotide in the reconstructed sequence |
 | `coverage_ratio` | The ratio of the reference sequence covered by the sequence of interest. |
 | `drug_class` | Set of antibiotic molecules, with similar chemical structures, molecular targets, and/or modes and mechanisms of action |
 | `input_gene_length` | The length (number of positions) of a target gene sequence submitted by a user |
 | `input_gene_start` | The position of the first nucleotide in a gene sequence being analyzed (input gene sequence) |
 | `input_gene_stop` | The position of the last nucleotide in a gene sequence being analyzed (input gene sequence) |
 | `input_protein_length` | The length (number of positions) of a protein target sequence submitted by a user |
 | `input_protein_start` | The position of the first amino acid in a protein sequence being analyzed (input protein sequence) |
 | `input_protein_stop` | The position of the last amino acid in a protein sequence being analyzed (input protein sequence) |
 | `input_sequence_id` | An identifier of molecular sequence(s) or entries from a molecular sequence database |
 | `nucleotide_mutation` | The nucleotide sequence change(s) detected in the sequence being analyzed compared to a reference |
 | `nucleotide_mutation_interpretation` | A description of identified nucleotide mutation(s) that facilitate clinical interpretation |
 | `predicted_phenotype` | A characteristic of an organism that is predicted rather than directly measured or observed |
 | `predicted_phenotype_confidence_level` | The level of confidence in a predicted phenotype |
 | `amino_acid_mutation` | The amino acid sequence change(s) detected in the sequence being analyzed compared to a reference |
 | `amino_acid_mutation_interpretation` | A description of identified amino acid mutation(s) that facilitate clinical interpretation. |
 | `reference_gene_length` | The length (number of positions) of a gene reference sequence retrieved from a database |
 | `reference_gene_start` | The position of the first nucleotide in a reference gene sequence |
 | `reference_gene_stop` | The position of the last nucleotide in a reference sequence |
 | `reference_protein_length` | The length (number of positions) of a protein reference sequence retrieved from a database |
 | `reference_protein_start` | The position of the first amino acid in a reference protein sequence |
 | `reference_protein_stop` | The position of the last amino acid in a reference protein sequence |
 | `resistance_mechanism` | Antibiotic resistance mechanisms evolve naturally via natural selection through random mutation, but it could also be engineered by applying an evolutionary stress on a population |
 | `strand_orientation` | The orientation of a genomic element on the double-stranded molecule |
 | `sequence_identity` | Sequence identity is the number (%) of matches (identical characters) in positions from an alignment of two molecular sequences |

</details>

[hAMRonization](https://github.com/pha4ge/hAMRonization) summarizes the outputs of the **antimicrobial resistance gene** detection tools (ABRicate, AMRFinderPlus, DeepARG, fARGene, RGI) into a single unified format. It supports a variety of summary options including an interactive summary.


### BGC detection tools:

[antiSMASH](#antismash), [deepBGC](#deepbgc), [GECCO](#gecco), [hmmsearch](#hmmsearch), summary: [comBGC](#combgc)

#### antiSMASH

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

#### deepBGC

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

#### GECCO

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

#### comBGC

<details markdown="1">
<summary>Output files</summary>

- `comBGC/`
  - `combgc_complete_summary.csv`: summarized output from all BGC detection tools used in tsv format.

</details>

<details markdown="1">
<summary>BGC summary table headers</summary>

| Table column | Description |
|---|---|
| `Sample_ID` | IDs of your samples |
| `Prediction_tool` | BGC prediction tool (antiSMASH, DeepBGC, and/or GECCO) |
| `Contig_ID` | ID of the contig containing the candidate BGC |
| `Product_class` | Predicted BGC type/class |
| `BGC_probability`| Confidence of BGC candidate as inferred by the respective tool |
| `BGC_complete` | Whether BGC sequence is assumed to be complete or truncated by the edge of the contig |
| `BGC_start` | Predicted BGC start position on the ontig |
| `BGC_end` | Predicted BGC end position on the contig |
| `BGC_length` | Length of the predicted BGC |
| `CDS_ID` | ID of the coding sequence(s) (CDS) from the annotation step (prodigal/prokka/bakta) if provided by the tool |
| `CDS_count` | Number of CDSs the BGC contains |
| `PFAM_domains` | Inferred PFAM IDs or annotations if provided by the tool |
| `MIBiG_ID` | Inferred MIBiG IDs if provided by the tool |
| `InterPro_ID` | Inferred InterPro IDs if provided by the tool |

</details>

**comBGC** is a local module which summarizes the results of the **Biosynthetic Gene Cluster (BGC)** prediction tools (antiSMASH, deepBGC, GECCO) used in the pipeline into one comprehensive summary with standardized headers.

**Note:** comBGC does not feature hmmer_hmmsearch at the moment.

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
