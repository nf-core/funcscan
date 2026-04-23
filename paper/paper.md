---
title: "nf-core/funcscan: A Nextflow pipeline to identify the biosynthetic potential and resistome of bacterial (meta)genomes"
tags:
  - nf-core
  - nextflow
  - pipeline
  - bioinformatics
  - AMP
  - AMR
  - antibiotic-resistance
  - antimicrobial-peptides
  - antimicrobial-resistance-genes
  - ARG
  - assembly
  - BGC
  - biosynthetic-gene-clusters
  - contigs
  - function
  - metagenomics
  - natural-products
  - screening
  - secondary-metabolites
  - taxonomic-classification
  - carbohydrate-active-enzyme
  - CAZyme-gene-cluster
  - CGC
authors:
  - name: Jasmin Frangenberg
    orcid: 0009-0004-5961-4709
    affiliation: 1
  - name: James A. Fellows Yates
    orcid: 0000-0001-5585-6277
    affiliation: "1, 2, 3" # (Multiple affiliations must be quoted)
  - name: Anan Ibrahim
    orcid: 0000-0003-3719-901X
    affiliation: 1
  - name: Louisa Perelo
    orcid: 0009-0002-6815-8608
    affiliation: 4
  - name: Haidong Yi
  - orcid:
    affiliation:
  - name: Xinpeng Zhang
  - orcid:
    affiliation:
  - name: Alexandru Mizeranschi
  - orcid:
    affiliation:
  - name: Dediu Codrin
  - orcid:
    affiliation:
  - name: Moritz E. Beber
    orcid: 0000-0003-2406-1978
    affiliation: 5
  - name: nf-core community
    affiliation: 6
  - name: Sven Nahnsen
    orcid: 0000-0002-4375-0691
    affiliation: "4, 7, 8"
  - name: Pierre Stallforth
    orcid: 0000-0001-7260-9921
    affiliation: "1, 11"
  - name: Christina Warinner
    orcid: 0000-0002-4528-5877
    affiliation: "2, 3, 9, 10"
affiliations:
  - name: Department of Paleobiotechnology, Leibniz Institute for Natural Product Research and Infection Biology Hans Knöll Institute, Germany
    index: 1
  - name: Department of Archaeogenetics, Max Planck Institute for Evolutionary Anthropology, Germany
    index: 2
  - name: Associated Research Group of Archaeogenetics, Leibniz Institute for Natural Product Research and Infection Biology Hans Knöll Institute, Germany
    index: 3
  - name: Quantitative Biology Center (QBiC), University of Tübingen, Germany
    index: 4
  - name: Institute for Globally Distributed Open Research and Education (IGDORE), Sweden
    index: 5
  - name: nf-core community members are available at acknowledgments.
    index: 6
  - name: M3 Research Center, Faculty of Medicine, University of Tübingen, Germany
    index: 7
  - name: Department of Computer Science, Institute for Bioinformatics and Medical Informatics (IBMI), University of Tübingen, Tübingen, Germany
    index: 8
  - name: Faculty of Biological Sciences, Friedrich-Schiller University Jena, Germany
    index: 9
  - name: Department of Anthropology, Harvard University, USA
    index: 10
  - name: Institute of Organic and Macromolecular Chemistry, Friedrich Schiller University Jena, Germany
    index: 11
date: 14 April 2026
bibliography: paper.bib
---

# Summary

Genome-mining of bacterial DNA fosters the discovery of antimicrobial resistance-related genes as well as genes required for the biosynthesis of low molecular weight natural products or specialised metabolites.
Despite the availability of many bioinformatic tools to identify such functional genes, screening of genomic features remains inefficient due to heterogeneous computational platforms, accessibility, scalability, and inconsistent reporting and formatting of the results.
Here, we present nf-core/funcscan, an open source bioinformatics pipeline for the screening of microbial functional features from assembled contigs or genomes.
The pipeline currently integrates 13 tools to simultaneously predict antimicrobial peptides, antibiotic resistance genes, biosynthetic gene clusters, and taxonomic classification from partial or full genomes.
It also introduces standardised and aggregated output file reports across all tools, enabling the rapid evaluation, visualisation, and interpretation of results.
Written in the Nextflow workflow language, it is straightforward to install, portable across platforms ranging from personal laptops to high-performance computing clusters, and fully reproducible via the use of software containers.

# Statement of need

The emergence and spread of multidrug resistant microbial pathogens poses a serious threat to global health.
Traditionally, most antiinfective drugs have been derived from bacterially produced low molecular weight natural products.
To ensure self-resistance against antimicrobial agents, the producing bacteria typically exhibit resistance mechanisms.
As a consequence, the evolution of antimicrobials and the corresponding resistance mechanisms are strongly correlated.
Although antibiotic resistance is tightly linked to self-protection of the producing organisms, the recent excessive use of antibiotics and lack of global surveillance both in healthcare and agriculture has led to an explosion of multidrug resistant bacteria.
Over the past few decades, the spread of antibiotic resistance genes (ARGs) and pathogenic bacteria carrying them has grown to a major threat to human health.
However, investigating antibiotic agents in combination with antibiotic resistance mechanisms and ARG evolution has the potential to aid in the development of new antibiotics.

Due to this pressing problem, a large suite of different tools has been developed for the rapid identification of different functional gene types.
These tools use different search algorithms and databases (e.g. deepBGC: machine-learning, antiSMASH: rule-based) for the prediction of microbial metabolites, which differ in quality and quantity of the predicted properties.
Thus, to maximise the potential of detecting important functional genes, researchers often need to use multiple approaches to ensure maximum detection sensitivity across all metabolite categories.
Since these tools are often developed as standalone tools they have to be executed separately.
This renders analyses inefficient and thus impedes scalability and poses the risk of lowering reproducibility.
While some tools are available as software containers (e.g. via docker, singularity), thus helping reproducibility of results, they require a series of steps to prepare input data and manually store and filter results.
Additionally, standalone tools have their own unique output format, which These points strongly hamper efficiency and in many cases reproducibility of complex analyses.
Overall, in order to obtain results from various tools in a uniform format, manual inspection is still necessary.
This renders the comparison of results from large datasets against multiple tools very impractical if not impossible.
Previous efforts to scale up the predictive power of different tools, including assembly, open reading frame (ORF) annotation, or gene-identifcation for functional prediction resulted in the genrationt of pipelines mettannotator, bacannot, SqueezeMeta, MetaErg, METABOLIC, HT-ARGfinder, ARGs-OAP, PathoFact, and antiSMASH.
However, so far, no pipeline has been created that allows for the identification and prediction of antimicrobial peptide (AMP) genes, ARGs, and biosynthetic gene clusters (BGCs) simultaneously from multiple samples in a harmonised manner.
Finally, extensive command-line knowledge, the use of shell scripts, and manual installation of software dependencies to run many of these tools, effectively precludes thir use by biochemists, biomolecular scientists, and biologists who typically have limited computational training.
Here, we present nf-core/funcscan, a Nextflow pipeline following nf-core best practices for the simultaneous screening of multiple functional and biosynthetic components from assembled contiguous sequences (contigs), specifically predicting ARGs, BGCs, AMP-encoding genes, and providing taxonomic information of the producing organisms from (meta)genomic sequences in a portable, reproducible, and scalable manner.
This allows researchers to obtain a holistic view on the genomic context of identified genes for downstream analyses in the context of antimicrobial resistance

# State of the field

The continuing decrease in sequencing costs and the subsequent increase in available sequenced prokaryotic genomes and metagenomes has gone hand-in-hand with the development of numerous bioinformatics tools to predict gene functions.
Several pipelines have been developed to chain single-purpose tools together to provide a more comprehensive context.
nf-core/funcscan is to the best of our knowledge the first pipeline to combine screening for AMPs, ARGs, and BGCs.
However, pipelines with similar functionality have been developed, with the closest one being mettannotator (Table 1).
This pipeline meets all criteria of scalability and reproducibility on the same level as nf-core/funcscan because it is likewise written in Nextflow and in most parts based on the nf-core pipeline template.
While focussing on somewhat different gene types (e.g. snRNA, mobilome), shared features include ARG and BGC prediction as well as aggregation of results.
In contrast, nf-core/funcscan provides additional AMP screening, CAZyme screening, and the integration of taxonomic classifications for all genes.
Regarding pipeline stability and reliability, nf-core/funcscan is the only pipeline to implement comprehensive unit tests on module and pipeline level, using the nf-test framework (Table \ref{tab:pipelines}).

| Feature                                 | funcscan | mettannotator | bacannot | HT-ARGfinder | PathoFact | SqueezeMeta | MetaERG | ARGs-OAP |
| --------------------------------------- | -------- | ------------- | -------- | ------------ | --------- | ----------- | ------- | -------- |
| ARG screening                           | ✓        | ✓             | ✓        | ✓            | ✓         | (✓)         | (✓)     | ✓        |
| AMP screening                           | ✓        | ✗             | ✗        | ✗            | ✗         | (✓)         | (✓)     | ✗        |
| BGC screening                           | ✓        | ✓             | ✗        | ✗            | ✗         | (✗)         | (✗)     | ✗        |
| CAZyme screening                        | ✓        | ✓             | ✗        | ✗            | ✗         | ✗           | ✗       | ✗        |
| Taxonomic assignment of contigs         | ✓        | ✗             | ✗        | (✗)          | (✗)       | ✓           | ✓       | ✗        |
| Results summary                         | ✓        | ✓             | ✓        | (✓)          | (✓)       | ✓           | ✓       | ✗        |
| Container support (docker, singularity) | ✓        | ✓             | ✓        | ✗            | ✗         | ✗           | ✓       | (✗)      |
| Modularity                              | ✓        | ✓             | ✓        | ✗            | ✓         | (✓)         | ✗       | ✗        |
| One-click installation                  | ✓        | ✓             | ✓        | ✗            | ✗         | (✗)         | ✗       | ✗        |
| Local installation possible             | ✓        | ✓             | ✓        | ✓            | ✓         | ✓           | ✓       | ✗        |
| Web-based execution possible            | (✓)      | (✓)           | (✓)      | ✗            | ✗         | ✗           | ✗       | ✗        |
| Software reviewing                      | ✓        | ✓             | ✗        | ✗            | ✗         | ✗           | ✗       | ✗        |
| Automated unit tests                    | ✓        | (✗)           | (✗)      | (✗)          | (✗)       | ✗           | ✗       | ✗        |
| License                                 | MIT      | Apache-2.0    | GPL-3.0  | None         | GPL-3.0   | GPL-3.0     | AFL     | AFL      |

: Comparison of nf-core/funcscan with other related pipelines for ARG, AMP, and BGC discovery. Parentheses indicate either unspecific gene screening or partly fulfilled criteria. \label{tab:pipelines}

# Software design

nf-core/funcscan simultaneously predicts antimicrobial peptide (AMP) genes, antibiotic resistance genes (ARGs), biosynthetic gene clusters (BGCs) as well as carbohydrate active enzyme gene clusters (CGC) from partial or full (meta)genomic sequences.
In addition, the bacterial taxonomy of input sequences is determined and standardised summaries of all tool outputs are provided (Fig. \ref{fig:workflow}).

![Workflow overview of nf-core/funcscan.
(1), genomic sequences are prepared and annotated with one of four ORF annotation tools.
Two additional classification workflows can be used to classify contigs taxonomically (light gray) or obtain additional protein domain information (dark gray).
(2), depending on which workflows are selected by the user, the biosynthetic gene cluster (BGC, purple), antimicrobial peptide (AMP, orange), antibiotic resistance gene (ARG, yellow), or carbohydrate-active enzymes (CAZyme) workflows with their customisable parameters are executed.
(3), the results of all tools for each gene category are aggregated and saved in a human- and machine-readable tabular format.\label{fig:workflow}](figure1.png)

## Input pre-processing and open reading frame annotation

The pipeline processes a two- to four- column table (comma-separated, CSV format) samplesheet as input.
Sample names and paths to the respective nucleotide FASTA files containing (meta)genomic contigs or genomes to be screened are required.
Optionally, preannotated sequence files can be supplied to the pipeline in a four-column samplesheet with open reading frame amino acid sequences in FASTA format, and their respective annotations in GenBank Flat File format.
During preprocessing, any gzipped sequence files are decompressed, and, when running the BGC subworkflow, short contigs are removed by SeqKit (default: contigs shorter than 3,000 bp).
The latter step reduces the runtime of the pipeline by removing too-short sequences that would produce no biologically meaningful BGC results.
Open reading frames are predicted from the pre-processed sequences by one of four annotation tools (Bakta, Prodigal, Prokka, and Pyrodigal).
If annotated sequence files as described above are provided in the samplesheet, this step is skipped.

Various tools of nf-core/funcscan rely on databases and reference files to operate.
The pipeline offers the functionality to download these databases automatically for the user, which can then be stored and reused in future pipeline runs to minimise pipeline runtime, network traffic, and possible download limits.
The database download is applicable for MMSeqs2, Bakta, AMPcombi, AMRFinderPlus, DeepARG, RGI, antiSMASH, DeepBGC, and InterProScan.

## Gene prediction and taxonomic classification

In a second step, users can choose to scan genomic sequences in parallel with three dedicated workflows for AMPs, ARGs and BGCs, applying up to currently a total of 12 gene identification tools:

- ARG subworkflow: ABRicate, AMRFinderPlus, DeepARG, fARGene, RGI
- BGC subworkflow: antiSMASH, DeepBGC, GECCO, hmmsearch
- AMP subworkflow: ampir, AMPlify, hmmsearch, Macrel

In an additional optional parallel screening step, all input sequences can be taxonomically classified by MMSeqs2 to determine likely source hosts of each functional hit.
Characterising the taxonomic origin of metagenomic contigs can inform users about potentially suitable hosts for downstream experiments, e.g. heterologous expression systems.
The taxonomic classification supports a variety of reference databases (e.g. GTDB, UniProt, UniRef, NR, Kalamari) to suit different user requirements.
Optionally, protein domains and families can be further annotated by InterProScan(37, 38).

## Aggregation of screening results

All screening tools of nf-core/funcscan have heterogeneous output formats and label their respective gene predictions differently.
This hampers aggregation and cross-comparisons of results, requiring manual inspection and ‘clean up’ of results for downstream interpretation.
To enable users to easily extract information for downstream analyses, nf-core/funcscan aggregates the output of all gene and taxonomic screening tools in each executed subworkflow into single human- and machine-readable tables in CSV format per gene type, thereby allowing direct comparison of gene classification with possible taxonomic sources.
For the summary of ARGs we have used the existing hAMRonization software. Since similar tools do not exist for AMPs and BGCs, we developed two novel tools for the aggregation of these gene types.
comBGC and AMPcombi parse the results of BGC and AMP prediction tools and summarise them into single tables, respectively.
Furthermore, AMPcombi aligns the AMP hits against a reference AMP database for deeper functional classification.
To assist researchers in their choice of genes for testing in wet-lab heterologous expression systems, AMPcombi provides the ability to reduce false positive hits by additional post-screening filtering steps of AMP results.
Reasonable default parameters are set by the pipeline and can be adjusted by the user.
Finally, three local pipeline modules merge the gene summaries with taxonomy results.

## Reproducibility and scalability

All nf-core pipelines utilise software environments (conda) or containers (Docker, Singularity), which have the advantage of isolating the dependencies of all workflows from each other and rendering pipeline execution highly reproducible, portable, and platform-independent.
Each tool of nf-core/funcscan is automatically pulled from the respective container registry when executing a pipeline run.
The pipeline itself is easy to install as it has only few minimum dependencies (Nextflow itself, and one of Docker, Singularity, Podman, Shifter, Charliecloud, and conda).
The configuration of the pipeline to the underlying computing system requires knowledge of its software environment and hardware resources.
To facilitate easy configuration, nf-core provides already centralised configurations for more than 150 HPCs via the central nf-core/configs repository (https://github.com/nf-core/configs).
The performance of each pipeline run (including software versions of all applied tools, memory and CPU usage) is summarsed in HTML reports for all steps of all subworkflows for users to estimate future runtime and/or computational resources.

# Research impact statement

nf-core/funcscan has developed an active user community of scientific users and developers who continuously contribute ideas, bug reports and code via issues and pull requests on GitHub.
An exemplary case in point is the contribution of a whole new workflow (CAZyme screening) by new community members.
Discussions of pipeline as well as research domain related topics happen on the open-to-join nf-core workspace on the Slack platform, illustrating the public interest and pro-active efforts from scientific users to use, maintain, and improve the pipeline functionalities.

# AI usage disclosure

No generative AI tools were used in the development of this software, the writing
of this manuscript, or the preparation of supporting materials.

# Acknowledgements

We thank Vedanth Ramji for adding argNorm to the ARG subworkflow.
A full list of nf-core community members is available at https://nf-co.re/contributors/.
We thank Martin Klapper and Rosa Herbst for helpful feedback on relevant BGC and AMP properties during comBGC and AMPcombi development.
J.F. received a fellowship from the International Leibniz Research School (under the head of the Jena School for Microbial Communication, JSMC).

This project was funded by grants from the Werner Siemens Foundation (Paleobiotechnology to C.W. and P.S.) and the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation, under Germany’s Excellence Strategy – EXC 2051 – Project-ID 390713860 to C.W. and P.S.).
J.A.F.Y was funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) – project number 460129525 (NFDI4Microbiota, FlexFund project EnterArchaeo).
This work was supported by the BMBF-funded de.NBI Cloud within the German Network for Bioinformatics Infrastructure (de.NBI) (031A532B, 031A533A, 031A533B, 031A534A, 031A535A, 031A537A, 031A537B, 031A537C, 031A537D, 031A538A).

# References
