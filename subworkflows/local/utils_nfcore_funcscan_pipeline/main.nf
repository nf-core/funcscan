//
// Subworkflow with functionality specific to the nf-core/pipeline pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromSamplesheet("input")
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map {
            validateInputSamplesheet(it)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ it.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
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

    def citation_text       = [
        preprocessing_text,
        annotation_text,
        params.run_amp_screening ? amp_text : "",
        params.run_arg_screening ? arg_text : "",
        params.run_bgc_screening ? bgc_text : "",
        postprocessing_text,
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def preprocessing_text = "<li>Li, H. (2023). bioawk: BWK awk modified for biological data. Github. Retrieved July 12, 2023, from <a href=\"https://github.com/lh3/bioawk\">https://github.com/lh3/bioawk</a></li>"

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

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
