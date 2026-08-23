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
    affiliation: "1, 2"
  - name: James A. Fellows Yates
    orcid: 0000-0001-5585-6277
    affiliation: "1, 3, 4"
  - name: Anan Ibrahim
    orcid: 0000-0003-3719-901X
    affiliation: 1
  - name: Louisa Perelo
    orcid: 0009-0002-6815-8608
    affiliation: 5
  - name: Haidong Yi
    orcid: 0000-0002-2591-1922
    affiliation: 6
  - name: Xinpeng Zhang
    orcid: 0000-0001-7567-8973
    affiliation: 7
  - name: Octavian Codrin Dediu
    orcid: 0009-0006-9204-8870
    affiliation: 8
  - name: Alexandru Mizeranschi
    orcid: 0000-0002-1168-6285
    affiliation: "9,10"
  - name: Moritz E. Beber
    orcid: 0000-0003-2406-1978
    affiliation: 11
  - name: nf-core community
    affiliation: 12
  - name: Sven Nahnsen
    orcid: 0000-0002-4375-0691
    affiliation: "5, 13, 14"
  - name: Pierre Stallforth
    orcid: 0000-0001-7260-9921
    affiliation: "1, 15"
  - name: Christina Warinner
    orcid: 0000-0002-4528-5877
    affiliation: "3, 4, 16, 17"
affiliations:
  - name: Department of Paleobiotechnology, Leibniz Institute for Natural Product Research and Infection Biology Hans Knöll Institute, Germany
    index: 1
  - name: TIB – Leibniz Information Centre for Science and Technology, Germany
    index: 2
  - name: Department of Archaeogenetics, Max Planck Institute for Evolutionary Anthropology, Germany
    index: 3
  - name: Associated Research Group of Archaeogenetics, Leibniz Institute for Natural Product Research and Infection Biology Hans Knöll Institute, Germany
    index: 4
  - name: Quantitative Biology Center (QBiC), University of Tübingen, Germany
    index: 5
  - name: St. Jude Children's Research Hospital, USA
    index: 6
  - name: Nebraska Food for Health Center, Department of Food Science and Technology, University of Nebraska, USA
    index: 7
  - name: Faculty of Computer Science, West University of Timisoara, Romania
    index: 8
  - name: Research and Development Station for Bovine – Arad, Romania
    index: 9
  - name: Institute for Advanced Environmental Research, West University of Timisoara, Romania
    index: 10
  - name: Institute for Globally Distributed Open Research and Education (IGDORE), Sweden
    index: 11
  - name: nf-core community members are available at acknowledgments.
    index: 12
  - name: M3 Research Center, Faculty of Medicine, University of Tübingen, Germany
    index: 13
  - name: Department of Computer Science, Institute for Bioinformatics and Medical Informatics (IBMI), University of Tübingen, Germany
    index: 14
  - name: Institute of Organic and Macromolecular Chemistry, Friedrich Schiller University Jena, Germany
    index: 15
  - name: Faculty of Biological Sciences, Friedrich-Schiller University Jena, Germany
    index: 16
  - name: Department of Anthropology, Harvard University, USA
    index: 17
date: 14 April 2026
bibliography: paper.bib
header-includes:
  - \usepackage{rotating}
  - \usepackage{booktabs}
---

# Summary

Genome-mining of bacterial DNA enables the discovery of antimicrobial resistance-related genes, genes required for the biosynthesis of low molecular weight natural products, and other specialised metabolites.
However, execution of the multiple bioinformatic tools used in screening analyses remains inefficient due to heterogenous software interfaces, reporting, and formatting of the output files of similar tools, which limits scalability of such analyses.

nf-core/funcscan is a portable and reproducible open source Nextflow bioinformatics pipeline for the screening of microbial functional features from assembled contigs or genomes.
The pipeline executes up to 13 tools to simultaneously identify antimicrobial peptides, antibiotic resistance genes, biosynthetic gene clusters, carbohydrate-activate enzymes, and performs taxonomic classification of partial or full genomes.
To facilitate efficient results comparison and evaluation, it supports cross-tool output file standardisation and aggregation.

# Statement of need

Researchers often use multiple tools to ensure maximum detection sensitivity during genomic screening for potential gene candidates, as each tool uses different search algorithms and microbial metabolite databases.
However, heterogenous installation, inputs, and execution interfaces of these stand-alone tools impedes scalability, and decreases reproducibility due to user-error when executed manually.
Additionally, each tool often has its own unique output formats, making cross-comparison of results between tools and databases non-trivial, and again requiring inefficient manual postprocessing and inspection.

This necessity for manual execution and postprocessing of heterogenous outputs impacts the discovery of new drugs.
For example, antibiotics are typically derived from naturally evolved, bacterially-produced, low molecular weight natural products, and the discovery rate of novel molecules has seen recent plateauing.
In combination with an explosion in the evolution of multidrug resistant bacteria [@ventola_antibiotic_2015; @perry_prehistory_2016; @rascovan_exploring_2016], and a lack of global surveillance both in healthcare and agriculture, this is contributing to a major threat to global health [@murray_global_2022; @world_health_organization_global_2022].
Therefore high-throughput and scalable approaches are needed to allow the rapid identification of metabolites from novel sources, as well as live monitoring of the spread of antibiotic resistance within microbial populations.

Here, we present nf-core/funcscan, a Nextflow [@di_tommaso_nextflow_2017] pipeline following nf-core best practices [@ewels_nf-core_2020;@Langer2025-th] for the automated and in-parallel screening of different functional gene groups with multiple tools and databases.
The pipeline currently supports detection of antimicrobial peptide (AMPs) genes, antimicrobial resistance genes (ARGs), biosynthetic gene clusters (BGCs), and carbohydrate-active enzyme gene clusters (CGCs).

# State of the field

Previous efforts to scale up the predictive power of different tools for functional gene prediction include mettannotator [@gurbich_mettannotator_2025], bacannot [@almeida_scalable_2023], PathoFact [@de_nies_pathofact_2021] SqueezeMeta [@tamames_squeezemeta_2019], MetaErg [@dong_integrated_2019], and ARGs-OAP [@yin_args-oap_2022] (Table \ref{tab:pipelines}).
However, to our knowledge, these are typically focused on singular gene categories or groups (e.g. antimicrobial resistance), aim to be 'end-to-end' pipelines including read preprocessing and assembly, or do not provide important contextual information about the potential hits (such as taxonomic information).

In particular, the main factors that distinguish nf-core/funcscan from the most similar pipeline, mettannotator, are: support for metagenomic assembly input (rather than just genomes); automated taxonomic classification of contigs; more tools for ARG screening; AMP screening; and confirmed executable on other infrastructure than HPCs.
Gene types which nf-core/funcscan does not screen for due to its focus on AMPs, ARGs, and BGCs but mettannotator does are CRISPR arrays, antiphage defense, non-coding RNA, and pseudogenes.

\begin{sidewaystable}
\centering
\caption{Comparison of nf-core/funcscan with other related pipelines for ARG, AMP, and BGC discovery. Parentheses indicate either unspecific gene screening or partly fulfilled criteria.}
\label{tab:pipelines}
\begin{tabular}{l|l|l|l|l|l|l|l}
\toprule
Feature & funcscan & mettannotator & bacannot & PathoFact & SqueezeMeta & MetaERG & ARGs-OAP \\
\hline
ARG screening & + & + & + & + & (+) & (+) & + \\
AMP screening & + & − & − & − & (+) & (+) & − \\
BGC screening & + & + & − & − & (−) & (−) & − \\
CAZyme screening & + & + & − & − & − & − & − \\
Taxonomic assignment of contigs & + & − & − & (−) & + & + & − \\
Results summary & + & + & + & (+) & + & + & − \\
Container support (Docker, Singularity) & + & + & + & − & − & + & (−) \\
Modularity & + & + & + & + & (+) & − & − \\
One-click installation & + & + & + & − & (−) & − & − \\
Local installation possible & + & + & + & + & + & + & − \\
Web-based execution possible & (+) & (+) & (+) & − & − & − & − \\
Software reviewing & + & + & − & − & − & − & − \\
Automated unit tests & + & + & (−) & (−) & − & − & − \\
License & MIT & Apache-2.0 & GPL-3.0 & GPL-3.0 & GPL-3.0 & AFL & AFL \\
\bottomrule
\end{tabular}
\end{sidewaystable}

Extensive command-line knowledge and manual installation of software dependencies are also often required to run many of the existing pipelines.
This can preclude use by biochemists, biologists, etc. who typically have limited computational training.
In contrast, nf-core/funcscan aims to reduce complexity by screening from already assembled sequences, and end on aggregation of the screening results, through multiple methods for execution (Table \ref{tab:pipelines}).

# Workflow overview

nf-core/funcscan simultaneously predicts AMPs, ARGs, BGCs as well as CGCs from input partial or full (meta)genomic sequences.
Output files from the tools of each of the categories are aggregated and standardised for easy cross-comparison (Fig. \ref{fig:workflow}).

![Workflow overview of nf-core/funcscan.
(1), genomic sequences are prepared and annotated with one of four open reading frame annotation tools.
Two additional classification workflows can be used to classify contigs taxonomically (light gray) or obtain additional protein domain information (dark gray).
(2), depending on user-choice, the biosynthetic gene cluster (BGC, purple), antimicrobial peptide genes (AMP, orange), antibiotic resistance genes (ARG, yellow), or carbohydrate-active enzymes (CAZyme, blue) workflows with their customisable parameters are executed.
(3), the results of all tools for each gene category are aggregated and saved in a human- and machine-readable tabular format.\label{fig:workflow}](figure1.png)

## Input preprocessing and open reading frame annotation

The pipeline takes a two- to five-column tabular sample-sheet as input (comma-separated, CSV format).
This sample-sheet contains sample names, paths to (meta)genomic FASTA files and optionally pre-generated amino-acid FASTA, GFF, or GBK format annotation files.

Preprocessing steps reduce runtime by removing too-short sequences with Seqkit [@Shen2024-mg], when they may produce no biologically meaningful results.
Open reading frames are optionally predicted from the preprocessed sequences by one of four prokaryotic annotation tools: Bakta [@schwengers_bakta_2021], Prodigal [@Hyatt2010-yv], Prokka [@Seemann2014-ee], and Pyrodigal [@Larralde2022-uu].

When required, the pipeline downloads required screening-tool databases automatically for the user, and makes them available for future pipeline runs to minimise runtime and network traffic.

## Gene screening and taxonomic classification

Users choose to scan genomic sequences in parallel with up-to four dedicated subworkflows for AMPs, ARGs, BGCs, and CAZymes.
Up to a total of 13 gene identification tools can be applied:

- **ARGs**: ABRicate [@torsten_seemann_abricate_2020], AMRFinderPlus [@feldgarden_amrfinderplus_2021;@feldgarden_validating_2019], DeepARG [@arango-argoty_deeparg_2018], fARGene [@berglund_identification_2019], RGI [@alcock_card_2023]
- **BGCs**: antiSMASH [@blin_antismash_2025], DeepBGC [@hannigan_deep_2019], GECCO [@carroll_accurate_2021], hmmsearch [@eddy_accelerated_2011]
- **AMPs**: ampir [@fingerhut_ampir_2021], AMPlify [@li_models_2023, @li_amplify_2022], hmmsearch, Macrel [@santos-junior_macrel_2020]
- **CAZymes**: run_dbCAN [@zheng_dbcan3_2023]

To provide users information about potentially suitable hosts for downstream experiments, e.g. heterologous expression systems [@porse_biochemical_2018], an additional optional parallel workflow can taxonomically classify input contigs with MMSeqs2 [@mirdita_fast_2021].
Optionally, generic protein domains and families can be further annotated with InterProScan [@jones_interproscan_2014].

Pipeline parameters can be adjusted by userwritten- or nf-core GUI ([https://nf-co.re/launch](https://nf-co.re/launch))-generated Nextflow parameter files, or command-line arguments.

## Aggregation of screening results

nf-core/funcscan integrates dedicated tools to aggregate and standardise heterogenous output formats of multiple screening tools into a single human- and machine-readable tables in CSV format per gene type.
nf-core uses hAMRonization [@mendes_hamronization_2024] for ARGs, AMPcombi [@herbst_actifensin_2025] for AMPs, and a custom script 'comBGC' for BGC tool output aggregation.
These summaries are finally optionally complemented with results from the taxonomic classification workflow.

Building on the aggragation of screening results, two optional downstream analyses can be executed for the ARG and BGC workflows.
First, the ARG summary provided by hAMRonization can be further normalised and mapped to the antibiotic resistance ontology (ARO) by argNorm [@ugarcina_perovic_argnorm_2025].
This enhances ARG annotation by categorising drugs that ARGs confer resistance to.
Secondly, BGCs predicted by antiSMASH and GECCO can be clustered into Gene Cluster Families (GCFs) by BiG-SLiCE [@kautsar_big-slice_2021] to enable comparative analysis of biosynthetic diversity across samples.

## Reproducibility and scalability

All nf-core pipelines utilise software environments [from the Bioconda project, @Gruning2018-vr] or containers [e.g. Docker, Singularity, primarily from the Biocontainers project, @Da_Veiga_Leprevost2017-gl] for each integrated tool.
This provides the advantage of isolating the dependencies of all workflows from each other, thereby reducing installation problems.
The pipeline is thus easy to install with few minimum dependencies - Nextflow itself, and one of Nextflow-supported container/software environment management systems.
For further portability, nf-core provides integrated configurations for more than 150 institutional computational infrastructure (e.g. HPCs) via nf-core/configs ([https://nf-co.re/configs](https://nf-co.re/configs)).
Users on these infrastructure thus can run the pipelines with no-set up via a single parameter.

# Research impact statement

nf-core/funcscan has an active user community of scientific users and developers on the nf-core Slack and GitHub ([https://nf-co.re/join](https://nf-co.re/join)).
For example, the pipeline received the contribution of the CAZyme screening from community members outside of the original developers.
User discussions and support on the pipeline and on related research topics occur on the nf-core Slack workspace.
This illustrates the public interest and proactive efforts from scientific users to use, maintain, and improve the pipeline.
The pipeline is also already actively being used in research [e.g., @Tighe2024-fq, @Janak2026-ek, @Istanbullugil2026-fg, @Liepa2026-sw].

# AI usage disclosure

No generative AI tools were used in the development of this software, the writing of this manuscript, or the preparation of supporting materials.

# Acknowledgements

We thank Vedanth Ramji for adding argNorm to the ARG subworkflow.
A full list of nf-core community members is available at [https://nf-co.re/contributors/](https://nf-co.re/contributors/).
We thank Martin Klapper and Rosa Herbst for helpful feedback on relevant BGC and AMP properties during comBGC and AMPcombi development.
J.F. received a fellowship from the International Leibniz Research School (under the head of the Jena School for Microbial Communication, JSMC).
This work was conducted while J.F. was affiliated with the Department of Paleobiotechnology, Leibniz Institute for Natural Product Research and Infection Biology Hans Knöll Institute, Germany.

This project was funded by grants from the Werner Siemens Foundation (Paleobiotechnology to C.W. and P.S.) and the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation, under Germany’s Excellence Strategy – EXC 2051 – Project-ID 390713860 to C.W. and P.S.).
J.A.F.Y and C.W. were funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) – project number 460129525 (NFDI4Microbiota, FlexFund project EnterArchaeo).
J.A.F.Y and C.W. were supported by the Max Planck Society.
O.C.D. and A.E.M. were supported with computing infrastructure access by the West University of Timisoara (hpc.uvt.ro), acquired via grant no. 240/2020, ID 911 POC/398/1/1, financed by the European structural funds and Romanian government funds, and the project "Romanian Hub for Artificial Intelligence - HRIA", in the program PCID-IF/709/PCIDIF_P4/OP1/RSO1.6/PCIDIF_A12 - Action 4.1 - MySMIS 351416.
This work was supported by the de.NBI Cloud within the German Network for Bioinformatics Infrastructure (de.NBI) and ELIXIR-DE (Forschungszentrum Jülich and W-de.NBI-001, W-de.NBI-004, W-de.NBI-008, W-de.NBI-010, W-de.NBI-013, W-de.NBI-014, W-de.NBI-016, W-de.NBI-022).

# References
