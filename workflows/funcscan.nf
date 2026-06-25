/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                         } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                } from 'plugin/nf-schema'
include { paramsSummaryMultiqc            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText          } from '../subworkflows/local/utils_nfcore_funcscan_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { ANNOTATION                      } from '../subworkflows/local/annotation'
include { PROTEIN_ANNOTATION              } from '../subworkflows/local/protein_annotation'
include { AMP                             } from '../subworkflows/local/amp'
include { ARG                             } from '../subworkflows/local/arg'
include { BGC                             } from '../subworkflows/local/bgc'
include { CAZYME                          } from '../subworkflows/local/cazyme'
include { TAXA_CLASS                      } from '../subworkflows/local/taxa_class'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { GUNZIP as GUNZIP_INPUT_PREP     } from '../modules/nf-core/gunzip'
include { SEQKIT_SEQ as SEQKIT_SEQ_LENGTH } from '../modules/nf-core/seqkit/seq'
include { SEQKIT_SEQ as SEQKIT_SEQ_FILTER } from '../modules/nf-core/seqkit/seq'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FUNCSCAN {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CONFIG FILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Some tools require uncompressed input
    ch_input_prep = ch_samplesheet
        .map { meta, fasta, faa, gbk, gff -> [meta + [category: 'all'], [fasta, faa, gbk, gff]] }
        .transpose()
        .branch {
            compressed: it[1].toString().endsWith('.gz')
            uncompressed: it[1]
        }

    GUNZIP_INPUT_PREP(ch_input_prep.compressed)
    ch_versions = ch_versions.mix(GUNZIP_INPUT_PREP.out.versions)

    // Merge all the already uncompressed and newly compressed FASTAs here into
    // a single input channel for downstream
    ch_intermediate_input = GUNZIP_INPUT_PREP.out.gunzip
        .mix(ch_input_prep.uncompressed)
        .groupTuple()
        .map { meta, files ->
            def fasta_found = files.find { it.toString().tokenize('.').last().matches('fasta|fas|fna|fa') }
            def faa_found = files.find { it.toString().endsWith('.faa') }
            def gff_found = files.find { it.toString().tokenize('.').last().matches('gff|gff3') }
            def gbk_found = files.find { it.toString().tokenize('.').last().matches('gbk|gbff') }
            def fasta = fasta_found != null ? fasta_found : []
            def faa = faa_found != null ? faa_found : []
            def gff = gff_found != null ? gff_found : []
            def gbk = gbk_found != null ? gbk_found : []

            [meta, fasta, faa, gff, gbk]
        }
        .branch { meta, fasta, faa, gff, gbk ->
            preannotated: gff != [] || gbk != []
            fastas: true
        }

    // Duplicate and filter the duplicated file for long contigs only for BGC
    // This is to speed up BGC run and prevent 'no hits found'  fails
    if (params.run_bgc_screening) {
        SEQKIT_SEQ_LENGTH(ch_intermediate_input.fastas.map { meta, fasta, faa, gff, gbk -> [meta, fasta] })
        ch_input_for_annotation = ch_intermediate_input.fastas
            .map { meta, fasta, protein, gff, gbk -> [meta, fasta] }
            .mix(SEQKIT_SEQ_LENGTH.out.fastx.map { meta, fasta -> [meta + [category: 'long'], fasta] })
            .filter { meta, fasta ->
                if (fasta != [] && fasta.isEmpty()) {
                    log.warn("[nf-core/funcscan] Sample ${meta.id} does not have contigs longer than ${params.bgc_mincontiglength} bp. Will not be screened for BGCs.")
                }
                !fasta.isEmpty()
            }
        ch_versions = ch_versions.mix(SEQKIT_SEQ_LENGTH.out.versions)
    }
    else {
        ch_input_for_annotation = ch_intermediate_input.fastas.map { meta, fasta, protein, gff, gbk -> [meta, fasta] }
    }

    /*
        ANNOTATION
    */

    // Some tools require annotated FASTAs
    if ((params.run_arg_screening && !params.arg_skip_deeparg) || params.run_amp_screening || params.run_bgc_screening || params.run_cazyme_screening) {
        ANNOTATION(ch_input_for_annotation)
        ch_versions = ch_versions.mix(ANNOTATION.out.versions)

        ch_new_annotation = ch_input_for_annotation
            .join(ANNOTATION.out.faa)
            .join(ANNOTATION.out.gff)
            .join(ANNOTATION.out.gbk)
    }
    else {
        ch_new_annotation = ch_intermediate_input.fastas
    }

    // Mix back the preannotated samples with the newly annotated ones
    ch_new_annotation_short = ch_new_annotation.filter { meta, fasta, faa, gff, gbk -> meta.category != 'long' }

    // Add gff_type to meta for cazyme screening
    if ((params.run_cazyme_screening && !params.cazyme_skip_dbcan && (!params.dbcan_skip_cgc || !params.dbcan_skip_substrate)) && params.annotation_tool in ['pyrodigal', 'prodigal', 'prokka', 'bakta']) {
        ch_new_annotation_for_mixing = ch_new_annotation_short.map { meta, fasta, faa, gff, gbk ->
            def new_meta = meta + [gff_type: 'prodigal']
            // Only Use 'prodigal' as dbcan does not distinguish 'pyrodigal' and 'prodigal'
            [new_meta, fasta, faa, gff, gbk]
        }
    }
    else {
        ch_new_annotation_for_mixing = ch_new_annotation_short
    }

    ch_prepped_input = ch_new_annotation_for_mixing
        .mix(ch_intermediate_input.preannotated)
        .multiMap { meta, fasta, faa, gff, gbk ->
            fastas: [meta, fasta]
            faas: [meta, faa]
            gffs: [meta, gff]
            gbks: [meta, gbk]
        }

    if (params.run_bgc_screening) {

        ch_prepped_input_long = ch_new_annotation
            .filter { meta, fasta, faa, gff, gbk -> meta.category == 'long' }
            .mix(ch_intermediate_input.preannotated)
            .multiMap { meta, fasta, faa, gff, gbk ->
                fastas: [meta, fasta]
                faas: [meta, faa]
                gffs: [meta, gff]
                gbks: [meta, gbk]
            }
    }

    /*
        TAXONOMIC CLASSIFICATION
    */

    // The final subworkflow reports need taxonomic classification.
    // This can be either on NT or AA level depending on annotation.
    // TODO: Only NT at the moment. AA tax. classification will be added only when its PR is merged.
    if (params.run_taxa_classification) {
        TAXA_CLASS(ch_prepped_input.fastas)
        ch_versions = ch_versions.mix(TAXA_CLASS.out.versions)
        ch_taxonomy_tsv = TAXA_CLASS.out.sample_taxonomy
    }
    else {

        ch_mmseqs_db = channel.empty()
        ch_taxonomy_querydb = channel.empty()
        ch_taxonomy_querydb_taxdb = channel.empty()
        ch_taxonomy_tsv = channel.empty()
    }

    /*
        PROTEIN ANNOTATION
    */
    if (params.run_protein_annotation) {
        def filtered_faas = ch_prepped_input.faas.filter { meta, file ->
            if (file != [] && file.isEmpty()) {
                log.warn("[nf-core/funcscan] Annotation of the following sample produced an empty FAA file. InterProScan classification of the CDS requiring this file will not be executed: ${meta.id}")
            }
            !file.isEmpty()
        }

        SEQKIT_SEQ_FILTER(filtered_faas)
        ch_versions = ch_versions.mix(SEQKIT_SEQ_FILTER.out.versions)
        ch_input_for_protein_annotation = SEQKIT_SEQ_FILTER.out.fastx

        PROTEIN_ANNOTATION(ch_input_for_protein_annotation)
        ch_versions = ch_versions.mix(PROTEIN_ANNOTATION.out.versions)

        ch_interproscan_tsv = PROTEIN_ANNOTATION.out.tsv.map { meta, file ->
            if (file == [] || file.isEmpty()) {
                log.warn("[nf-core/funcscan] Protein annotation with InterProScan produced an empty TSV file. No protein annotation will be added for sample ${meta.id}.")
                [meta, []]
            }
            else {
                [meta, file]
            }
        }
    }
    else {
        ch_interproscan_tsv = ch_prepped_input.faas.map { meta, _files ->
            [meta, []]
        }
    }


    /*
        SCREENING
    */

    /*
        AMPs
    */
    if (params.run_amp_screening && !params.run_taxa_classification) {
        AMP(
            ch_prepped_input.fastas,
            ch_prepped_input.faas.filter { meta, file ->
                if (file != [] && file.isEmpty()) {
                    log.warn("[nf-core/funcscan] Annotation of following sample produced an empty FAA file. AMP screening tools requiring this file will not be executed: ${meta.id}")
                }
                !file.isEmpty()
            },
            ch_taxonomy_tsv,
            ch_prepped_input.gbks,
            ch_interproscan_tsv,
        )
        ch_versions = ch_versions.mix(AMP.out.versions)
    }
    else if (params.run_amp_screening && params.run_taxa_classification) {
        AMP(
            ch_prepped_input.fastas,
            ch_prepped_input.faas.filter { meta, file ->
                if (file != [] && file.isEmpty()) {
                    log.warn("[nf-core/funcscan] Annotation of following sample produced an empty FAA file. AMP screening tools requiring this file will not be executed: ${meta.id}")
                }
                !file.isEmpty()
            },
            ch_taxonomy_tsv.filter { meta, file ->
                if (file != [] && file.isEmpty()) {
                    log.warn("[nf-core/funcscan] Taxonomy classification of the following sample produced an empty TSV file. Taxonomy merging will not be executed: ${meta.id}")
                }
                !file.isEmpty()
            },
            ch_prepped_input.gbks,
            ch_interproscan_tsv,
        )
        ch_versions = ch_versions.mix(AMP.out.versions)
    }

    /*
        ARGs
    */
    if (params.run_arg_screening && !params.run_taxa_classification) {
        if (params.arg_skip_deeparg) {
            ARG(
                ch_prepped_input.fastas,
                [],
                ch_taxonomy_tsv,
            )
        }
        else {
            ARG(
                ch_prepped_input.fastas,
                ch_prepped_input.faas.filter { meta, file ->
                    if (file.isEmpty()) {
                        log.warn("[nf-core/funcscan] Annotation of following sample produced an empty FAA file. ARG screening tools requiring this file will not be executed: ${meta.id}")
                    }
                    !file.isEmpty()
                },
                ch_taxonomy_tsv,
            )
        }
        ch_versions = ch_versions.mix(ARG.out.versions)
    }
    else if (params.run_arg_screening && params.run_taxa_classification) {
        if (params.arg_skip_deeparg) {
            ARG(
                ch_prepped_input.fastas,
                [],
                ch_taxonomy_tsv.filter { meta, file ->
                    if (file.isEmpty()) {
                        log.warn("[nf-core/funcscan] Taxonomy classification of the following sample produced an empty TSV file. Taxonomy merging will not be executed: ${meta.id}")
                    }
                    !file.isEmpty()
                },
            )
        }
        else {
            ARG(
                ch_prepped_input.fastas,
                ch_prepped_input.faas.filter { meta, file ->
                    if (file.isEmpty()) {
                        log.warn("[nf-core/funcscan] Annotation of following sample produced an empty FAA file. ARG screening tools requiring this file will not be executed: ${meta.id}")
                    }
                    !file.isEmpty()
                },
                ch_taxonomy_tsv.filter { meta, file ->
                    if (file.isEmpty()) {
                        log.warn("[nf-core/funcscan] Taxonomy classification of the following sample produced an empty TSV file. Taxonomy merging will not be executed: ${meta.id}")
                    }
                    !file.isEmpty()
                },
            )
        }
        ch_versions = ch_versions.mix(ARG.out.versions)
    }

    /*
        BGCs
    */
    if (params.run_bgc_screening && !params.run_taxa_classification) {
        BGC(
            ch_prepped_input_long.fastas,
            ch_prepped_input_long.faas.filter { meta, file ->
                if (file != [] && file.isEmpty()) {
                    log.warn("[nf-core/funcscan] Annotation of following sample produced an empty GFF file. BGC screening tools requiring this file will not be executed: ${meta.id}")
                }
                !file.isEmpty()
            },
            ch_prepped_input_long.gbks.filter { meta, file ->
                if (file != [] && file.isEmpty()) {
                    log.warn("[nf-core/funcscan] Annotation of following sample produced an empty FAA file. BGC screening tools requiring this file will not be executed: ${meta.id}")
                }
                !file.isEmpty()
            },
            ch_taxonomy_tsv,
        )
        ch_versions = ch_versions.mix(BGC.out.versions)
    }
    else if (params.run_bgc_screening && params.run_taxa_classification) {
        BGC(
            ch_prepped_input_long.fastas,
            ch_prepped_input_long.faas.filter { meta, file ->
                if (file.isEmpty()) {
                    log.warn("[nf-core/funcscan] Annotation of following sample produced an empty FAA file. BGC screening tools requiring this file will not be executed: ${meta.id}")
                }
                !file.isEmpty()
            },
            ch_prepped_input_long.gbks.filter { meta, file ->
                if (file.isEmpty()) {
                    log.warn("[nf-core/funcscan] Annotation of following sample produced an empty GBK file. BGC screening tools requiring this file will not be executed: ${meta.id}")
                }
                !file.isEmpty()
            },
            ch_taxonomy_tsv.filter { meta, file ->
                if (file.isEmpty()) {
                    log.warn("[nf-core/funcscan] Taxonomy classification of the following sample produced an empty TSV file. Taxonomy merging will not be executed: ${meta.id}")
                }
                !file.isEmpty()
            },
        )
        ch_versions = ch_versions.mix(BGC.out.versions)
    }

    /*
        CAZYMEs
    */
    if (params.run_cazyme_screening) {
        CAZYME(
            ch_prepped_input.faas.filter { meta, file ->
                if (file != [] && file.isEmpty()) {
                    log.warn("[nf-core/funcscan] Annotation of following sample produced an empty FAA file. CAZyme screening tools requiring this file will not be executed: ${meta.id}")
                }
                !file.isEmpty()
            },
            ch_prepped_input.gffs,
        )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_' + 'funcscan_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    if ((params.run_arg_screening && !params.arg_skip_deeparg) || (params.run_amp_screening && (params.amp_run_hmmsearch || !params.amp_skip_amplify || !params.amp_skip_ampir)) || params.run_bgc_screening) {
        ch_multiqc_files = ch_multiqc_files.mix(ANNOTATION.out.multiqc_files.collect { it[1] })
    }

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'funcscan'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
