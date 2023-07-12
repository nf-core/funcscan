//
// This file holds several functions specific to the workflow/funcscan.nf in the nf-core/funcscan pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowFuncscan {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        genomeExistsError(params, log)

        //if (!params.fasta) {
        //    Nextflow.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
        //}
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Generate methods description for MultiQC
    //

    public static String toolCitationText(params) {

        def preprocessing_text = "The pipeline used the following tools: preprocessing included bioawk (Li 2023)."

        def annotation_text    = [
                "Annotation was carried out with:",
                params.annotation_tool == 'prodigal'  ? "Prodigal (Hyatt et al. 2010)." : "",
                params.annotation_tool == 'pyrodigal' ? "Pyrodigal (Larralde 2022)." : "",
                params.annotation_tool == 'bakta'     ? "BAKTA (Schwengers et al. 2021)." : "",
                params.annotation_tool == 'prokka'    ? "PROKKA (Seemann 2014)." : "",
            ].join(' ').trim()

        def amp_text           = [
                "The following antimicrobial peptide screening tools were used:",
                !params.amp_skip_amplify ? "AMPlify (Li et al. 2022)," : "",
                !params.amp_skip_macrel ? "Macrel (Santos-Júnior et al. 2020)," : "",
                !params.amp_skip_ampir ? "ampir (Fingerhut et al. 2021)," : "",
                !params.amp_skip_hmmsearch ? "HMMER (Eddy 2011)," : "",
                ". The output from the antimicrobial peptide screening tools were standardised and summarised with AMPcombi (Ibrahim and Perelo 2023)."
            ].join(' ').trim().replaceAll(", \\.", ".")

        def arg_text           = [
                "The following antimicrobial resistance gene screening tools were used:",
                !params.arg_skip_fargene ? "fARGene (Berglund et al. 2019)," : "",
                !params.arg_skip_rgi ? "RGI (Alcock et al. 2020)," : "",
                !params.arg_skip_amrfinderplus ? "AMRfinderplus (Feldgarden et al. 2021)," : "",
                !params.arg_skip_deeparg ? "deepARG (Arango-Argoty 2018)," : "",
                !params.arg_skip_abricate ? "ABRicate (Seemann 2020)," : "",
                ". The output from the antimicrobial resistance gene screening tools were standardised and summarised with hAMRonization (Maguire et al. 2023)."
            ].join(' ').trim().replaceAll(", +\\.", ".")

        def bgc_text           = [
                "The following biosynthetic gene cluster screening tools were used:",
                !params.bgc_skip_antismash ? "antiSMASH (Blin et al. 2021)," : "",
                !params.bgc_skip_deepbgc ? "deepBGC (Hannigan et al. 2019)," : "",
                !params.bgc_skip_gecco ? "GECCO (Carroll et al. 2021)," : "",
                !params.bgc_skip_hmmsearch ? "HMMER (Eddy 2011)," : "",
                ". The output from the biosynthetic gene cluster screening tools were standardised and summarised with comBGC (Frangenberg et al. 2023)."
            ].join(' ').replaceAll(", +\\.", ".").trim()

        def postprocessing_text = "Run statistics were reported using MultiQC (Ewels et al. 2016)."

        def citation_text = [
            preprocessing_text,
            annotation_text,
            params.run_amp_screening ? amp_text : "",
            params.run_arg_screening ? arg_text : "",
            params.run_bgc_screening ? bgc_text : "",
            postprocessing_text,
        ].join(' ').trim()

        return citation_text
    }

    public static String toolBibliographyText(params) {

        def preprocessing_text = "<li>Li, H. (2023). bioawk: BWK awk modified for biological data. Github. Retrieved July 12, 2023, from <a href=\"https://github.com/lh3/bioawk\">https://github.com/lh3/bioawk</a></li>"


        // TODO DUPLCIATED HMMSEARCH
        def annotation_text    = [
                params.annotation_tool == 'prodigal'  ? "<li>Hyatt, D., Chen, G. L., Locascio, P. F., Land, M. L., Larimer, F. W., & Hauser, L. J. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC bioinformatics, 11, 119. DOI: <a href=\"https://doi.org/10.1186/1471-2105-11-119\">10.1186/1471-2105-11-119</a>" : "",
                params.annotation_tool == 'pyrodigal' ? "<li>Larralde, M. (2022). Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes. Journal of Open Source Software, 7(72), 4296. DOI: <a href=\"https://doi.org/10.21105/joss.04296\">10.21105/joss.04296</a></li>" : "",
                params.annotation_tool == 'bakta'     ? "<li>Schwengers, O., Jelonek, L., Dieckmann, M. A., Beyvers, S., Blom, J., & Goesmann, A. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. Microbial Genomics, 7(11). DOI: <a href=\"https://doi.org/10.1099/mgen.0.000685\">10.1099/mgen.0.000685</a></li>" : "",
                params.annotation_tool == 'prokka'    ? "<li>Seemann, T. (2014). Prokka: rapid prokaryotic genome annotation. Bioinformatics (Oxford, England), 30(14), 2068–2069. DOI: <a href=\"https://doi.org/10.1093/bioinformatics/btu153\">10.1093/bioinformatics/btu153</a></li>" : "",
            ].join(' ').trim()

        def amp_text           = [
                !params.amp_skip_amplify   ? "<li>Li, C., Sutherland, D., Hammond, S. A., Yang, C., Taho, F., Bergman, L., Houston, S., Warren, R. L., Wong, T., Hoang, L., Cameron, C. E., Helbing, C. C., & Birol, I. (2022). AMPlify: attentive deep learning model for discovery of novel antimicrobial peptides effective against WHO priority pathogens. BMC genomics, 23(1), 77. DOI: <a href=\"https://doi.org/10.1186/s12864-022-08310-4\">10.1186/s12864-022-08310-4</a></li>" : "",
                !params.amp_skip_macrel    ? "<li>Santos-Júnior, C. D., Pan, S., Zhao, X. M., & Coelho, L. P. (2020). Macrel: antimicrobial peptide screening in genomes and metagenomes. PeerJ, 8, e10555. DOI: <a href=\"https://doi.org/10.7717/peerj.10555\">10.7717/peerj.10555</a></li>" : "",
                !params.amp_skip_ampir     ? "<li>Fingerhut, L., Miller, D. J., Strugnell, J. M., Daly, N. L., & Cooke, I. R. (2021). ampir: an R package for fast genome-wide prediction of antimicrobial peptides. Bioinformatics (Oxford, England), 36(21), 5262–5263. DOI: <a href=\"https://doi.org/10.1093/bioinformatics/btaa653\">10.1093/bioinformatics/btaa653</a></li>" : "",
                "<li>Ibrahim, A. & Perelo, L. (2023). Darcy220606/AMPcombi. DOI: <a href=\"https://doi.org/10.5281/zenodo.7639121\">10.5281/zenodo.7639121</a></li>"
            ].join(' ').trim().replaceAll(", \\.", ".")

        def arg_text           = [
                !params.arg_skip_fargene       ? "<li>Berglund, F., Österlund, T., Boulund, F., Marathe, N. P., Larsson, D., & Kristiansson, E. (2019). Identification and reconstruction of novel antibiotic resistance genes from metagenomes. Microbiome, 7(1), 52. DOI: <a href=\"https://doi.org/10.1186/s40168-019-0670-1\">10.1186/s40168-019-0670-1</a></li>" : "",
                !params.arg_skip_rgi           ? "<li>Alcock, B. P., Raphenya, A. R., Lau, T., Tsang, K. K., Bouchard, M., Edalatmand, A., Huynh, W., Nguyen, A. V., Cheng, A. A., Liu, S., Min, S. Y., Miroshnichenko, A., Tran, H. K., Werfalli, R. E., Nasir, J. A., Oloni, M., Speicher, D. J., Florescu, A., Singh, B., Faltyn, M., … McArthur, A. G. (2020). CARD 2020: antibiotic resistome surveillance with the comprehensive antibiotic resistance database. Nucleic acids research, 48(D1), D517–D525. DOI: <a href=\"https://doi.org/10.1093/nar/gkz935\">10.1093/nar/gkz935</a></li>" : "",
                !params.arg_skip_amrfinderplus ? "<li>Feldgarden, M., Brover, V., Gonzalez-Escalona, N., Frye, J. G., Haendiges, J., Haft, D. H., Hoffmann, M., Pettengill, J. B., Prasad, A. B., Tillman, G. E., Tyson, G. H., & Klimke, W. (2021). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Scientific reports, 11(1), 12728. DOI: <a href=\"https://doi.org/10.1038/s41598-021-91456-0\">10.1038/s41598-021-91456-0</a></li>" : "",
                !params.arg_skip_deeparg       ? "<li>Arango-Argoty, G., Garner, E., Pruden, A., Heath, L. S., Vikesland, P., & Zhang, L. (2018). DeepARG: a deep learning approach for predicting antibiotic resistance genes from metagenomic data. Microbiome, 6(1), 23. DOI: <a href=\"https://doi.org/10.1186/s40168-018-0401-z\">10.1186/s40168-018-0401-z" : "",
                !params.arg_skip_abricate      ? "<li>Seemann, T. (2020). ABRicate. Github <a href=\"https://github.com/tseemann/abricate\">https://github.com/tseemann/abricate</a>.</li>" : "",
                "<li>Public Health Alliance for Genomic Epidemiology (pha4ge). (2022). Parse multiple Antimicrobial Resistance Analysis Reports into a common data structure. Github. Retrieved October 5, 2022, from <a href=\"https://github.com/pha4ge/hAMRonization\">https://github.com/pha4ge/hAMRonization</a></li>"
            ].join(' ').trim().replaceAll(", +\\.", ".")

        def bgc_text           = [
                !params.bgc_skip_antismash ? "<li>Blin, K., Shaw, S., Kloosterman, A. M., Charlop-Powers, Z., van Wezel, G. P., Medema, M. H., & Weber, T. (2021). antiSMASH 6.0: improving cluster detection and comparison capabilities. Nucleic acids research, 49(W1), W29–W35. DOI: <a href=\"https://doi.org/10.1093/nar/gkab335\"10.1093/nar/gkab335</a></li>" : "",
                !params.bgc_skip_deepbgc   ? "<li>Hannigan, G. D., Prihoda, D., Palicka, A., Soukup, J., Klempir, O., Rampula, L., Durcak, J., Wurst, M., Kotowski, J., Chang, D., Wang, R., Piizzi, G., Temesi, G., Hazuda, D. J., Woelk, C. H., & Bitton, D. A. (2019). A deep learning genome-mining strategy for biosynthetic gene cluster prediction. Nucleic acids research, 47(18), e110. DOI: <a href=\"https://doi.org/10.1093/nar/gkz654\">10.1093/nar/gkz654</a></li>" : "",
                !params.bgc_skip_gecco     ? "<li>Carroll, L. M. , Larralde, M., Fleck, J. S., Ponnudurai, R., Milanese, A., Cappio Barazzone, E. & Zeller, G. (2021). Accurate de novo identification of biosynthetic gene clusters with GECCO. bioRxiv DOI: <a href=\"https://doi.org/10.1101/2021.05.03.442509\">0.1101/2021.05.03.442509</a></li>" : "",
                "<li>Frangenberg, J. Fellows Yates, J. A., Ibrahim, A., Perelo, L., & Beber, M. E. (2023). nf-core/funcscan: 1.0.0 - German Rollmops - 2023-02-15. <a href=\"https://doi.org/10.5281/zenodo.7643100\">https://doi.org/10.5281/zenodo.7643100</a></li>"
            ].join(' ').replaceAll(", +\\.", ".").trim()

        def postprocessing_text = "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. <a href=\"https://doi.org/10.1093/bioinformatics/btw354\">https://doi.org/10.1093/bioinformatics/btw354</a></li>"

        // Special as reused in multiple subworkflows, and we don't want to cause duplicates
        def hmmsearch_text = ( params.run_amp_screening && !params.amp_skip_hmmsearch ) || (params.run_bgc_screening && !params.bgc_skip_hmmsearch) ? "<li>Eddy S. R. (2011). Accelerated Profile HMM Searches. PLoS computational biology, 7(10), e1002195. DOI: <a href=\"https://doi.org/10.1371/journal.pcbi.1002195\">10.1371/journal.pcbi.1002195</a></li>" : ""

        def reference_text = [
            preprocessing_text,
            annotation_text,
            params.run_amp_screening ? amp_text : "",
            params.run_arg_screening ? arg_text : "",
            params.run_bgc_screening ? bgc_text : "",
            hmmsearch_text,
            postprocessing_text,
        ].join(' ').trim()

        return reference_text

    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml, params) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        // Pipeline DOI
        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        // Tool references
        meta["tool_citations"] = ""
        meta["tool_bibliography"] = ""

        meta["tool_citations"] = toolCitationText(params).replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
        meta["tool_bibliography"] = toolBibliographyText(params)

        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }
}
