/*
    Run annotation tools
*/

include { PROKKA                         } from '../../modules/nf-core/prokka/main'
include { PRODIGAL                       } from '../../modules/nf-core/prodigal/main'
include { PYRODIGAL                      } from '../../modules/nf-core/pyrodigal/main'
include { BAKTA_BAKTADBDOWNLOAD          } from '../../modules/nf-core/bakta/baktadbdownload/main'
include { BAKTA_BAKTA                    } from '../../modules/nf-core/bakta/bakta/main'
include { GUNZIP as GUNZIP_PRODIGAL_FNA  } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PRODIGAL_FAA  } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PRODIGAL_GBK  } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PYRODIGAL_FNA } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PYRODIGAL_FAA } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PYRODIGAL_GBK } from '../../modules/nf-core/gunzip/main'

workflow ANNOTATION {
    take:
    fasta // tuple val(meta), path(contigs)

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if (params.annotation_tool == "pyrodigal" || (params.annotation_tool == "prodigal" && params.run_bgc_screening == true && (!params.bgc_skip_antismash || !params.bgc_skip_deepbgc || !params.bgc_skip_gecco)) || (params.annotation_tool == "prodigal" && params.run_amp_screening == true)) {
        // Need to use Pyrodigal for most BGC tools and AMPcombi because Prodigal GBK annotation format is incompatible with them.

        if (params.annotation_tool == "prodigal" && params.run_bgc_screening == true && (!params.bgc_skip_antismash || !params.bgc_skip_deepbgc || !params.bgc_skip_gecco)) {
            log.warn("[nf-core/funcscan] Switching annotation tool to: Pyrodigal. This is because Prodigal annotations (in GBK format) are incompatible with antiSMASH, DeepBGC, and GECCO. If you specifically wish to run Prodigal instead, please skip antiSMASH, DeepBGC, and GECCO or provide a pre-annotated GBK file in the samplesheet.")
        }
        else if (params.annotation_tool == "prodigal" && params.run_amp_screening == true) {
            log.warn("[nf-core/funcscan] Switching annotation tool to: Pyrodigal. This is because Prodigal annotations (in GBK format) are incompatible with AMPcombi. If you specifically wish to run Prodigal instead, please skip AMP workflow or provide a pre-annotated GBK file in the samplesheet.")
        }

        PYRODIGAL(fasta, "gbk")
        GUNZIP_PYRODIGAL_FAA(PYRODIGAL.out.faa)
        GUNZIP_PYRODIGAL_FNA(PYRODIGAL.out.fna)
        GUNZIP_PYRODIGAL_GBK(PYRODIGAL.out.annotations)
        ch_versions = ch_versions.mix(PYRODIGAL.out.versions)
        ch_versions = ch_versions.mix(GUNZIP_PYRODIGAL_FAA.out.versions)
        ch_versions = ch_versions.mix(GUNZIP_PYRODIGAL_FNA.out.versions)
        ch_versions = ch_versions.mix(GUNZIP_PYRODIGAL_GBK.out.versions)
        ch_annotation_faa = GUNZIP_PYRODIGAL_FAA.out.gunzip
        ch_annotation_fna = GUNZIP_PYRODIGAL_FNA.out.gunzip
        ch_annotation_gbk = GUNZIP_PYRODIGAL_GBK.out.gunzip
    }
    else if (params.annotation_tool == "prodigal") {

        PRODIGAL(fasta, "gbk")
        GUNZIP_PRODIGAL_FAA(PRODIGAL.out.amino_acid_fasta)
        GUNZIP_PRODIGAL_FNA(PRODIGAL.out.nucleotide_fasta)
        GUNZIP_PRODIGAL_GBK(PRODIGAL.out.gene_annotations)
        ch_versions = ch_versions.mix(PRODIGAL.out.versions)
        ch_versions = ch_versions.mix(GUNZIP_PRODIGAL_FAA.out.versions)
        ch_versions = ch_versions.mix(GUNZIP_PRODIGAL_FNA.out.versions)
        ch_versions = ch_versions.mix(GUNZIP_PRODIGAL_GBK.out.versions)
        ch_annotation_faa = GUNZIP_PRODIGAL_FAA.out.gunzip
        ch_annotation_fna = GUNZIP_PRODIGAL_FNA.out.gunzip
        ch_annotation_gbk = GUNZIP_PRODIGAL_GBK.out.gunzip
    }
    else if (params.annotation_tool == "prokka") {

        PROKKA(fasta, [], [])
        ch_versions = ch_versions.mix(PROKKA.out.versions)
        ch_multiqc_files = PROKKA.out.txt.collect { it[1] }.ifEmpty([])
        ch_annotation_faa = PROKKA.out.faa
        ch_annotation_fna = PROKKA.out.fna
        ch_annotation_gbk = PROKKA.out.gbk
    }
    else if (params.annotation_tool == "bakta") {

        // BAKTA prepare download
        if (params.annotation_bakta_db) {
            ch_bakta_db = Channel
                .fromPath(params.annotation_bakta_db, checkIfExists: true)
                .first()
        }
        else {
            BAKTA_BAKTADBDOWNLOAD()
            ch_versions = ch_versions.mix(BAKTA_BAKTADBDOWNLOAD.out.versions)
            ch_bakta_db = BAKTA_BAKTADBDOWNLOAD.out.db
        }

        BAKTA_BAKTA(fasta, ch_bakta_db, [], [])
        ch_versions = ch_versions.mix(BAKTA_BAKTA.out.versions)
        ch_multiqc_files = BAKTA_BAKTA.out.txt.collect { it[1] }.ifEmpty([])
        ch_annotation_faa = BAKTA_BAKTA.out.faa
        ch_annotation_fna = BAKTA_BAKTA.out.fna
        ch_annotation_gbk = BAKTA_BAKTA.out.gbff
    }

    emit:
    versions      = ch_versions
    multiqc_files = ch_multiqc_files
    faa           = ch_annotation_faa // [ [meta], path(faa) ]
    fna           = ch_annotation_fna // [ [meta], path(fna) ]
    gbk           = ch_annotation_gbk // [ [meta], path(gbk) ]
}
