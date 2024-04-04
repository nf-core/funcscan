/*
    Run annotation tools
*/

include { PROKKA                         } from '../../modules/nf-core/prokka/main'
include { PRODIGAL as PRODIGAL_GFF       } from '../../modules/nf-core/prodigal/main'
include { PRODIGAL as PRODIGAL_GBK       } from '../../modules/nf-core/prodigal/main'
include { PYRODIGAL                      } from '../../modules/nf-core/pyrodigal/main'
include { BAKTA_BAKTADBDOWNLOAD          } from '../../modules/nf-core/bakta/baktadbdownload/main'
include { BAKTA_BAKTA                    } from '../../modules/nf-core/bakta/bakta/main'
include { GUNZIP as GUNZIP_PRODIGAL_FNA  } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PRODIGAL_FAA  } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PRODIGAL_GFF  } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PRODIGAL_GBK  } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PYRODIGAL_FNA } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PYRODIGAL_FAA } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PYRODIGAL_GFF } from '../../modules/nf-core/gunzip/main'

workflow ANNOTATION {
    take:
    fasta // tuple val(meta), path(contigs)

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

        // For prodigal: run twice, once for gff and once for gbk generation, (for parity with PROKKA which produces both)
        if ( params.annotation_tool == "prodigal" ) {
            PRODIGAL_GFF ( fasta, "gff" )
            GUNZIP_PRODIGAL_FAA ( PRODIGAL_GFF.out.amino_acid_fasta )
            GUNZIP_PRODIGAL_FNA ( PRODIGAL_GFF.out.nucleotide_fasta)
            GUNZIP_PRODIGAL_GFF ( PRODIGAL_GFF.out.gene_annotations )
            ch_versions       = ch_versions.mix(PRODIGAL_GFF.out.versions)
            ch_annotation_faa = GUNZIP_PRODIGAL_FAA.out.gunzip
            ch_annotation_fna = GUNZIP_PRODIGAL_FNA.out.gunzip
            ch_annotation_gff = GUNZIP_PRODIGAL_GFF.out.gunzip
            ch_annotation_gbk = Channel.empty() // Prodigal GBK and GFF output are mutually exclusive

            if ( params.save_annotations == true ) {
                PRODIGAL_GBK ( fasta, "gbk" )
                ch_versions       = ch_versions.mix(PRODIGAL_GBK.out.versions)
                ch_annotation_gbk = PRODIGAL_GBK.out.gene_annotations // Prodigal GBK output stays zipped because it is currently not used by any downstream subworkflow.
            }

        } else if ( params.annotation_tool == "pyrodigal" ) {

            PYRODIGAL ( fasta, "gbk" )
            GUNZIP_PYRODIGAL_FAA ( PYRODIGAL.out.faa )
            GUNZIP_PYRODIGAL_FNA ( PYRODIGAL.out.fna)
            GUNZIP_PYRODIGAL_GFF ( PYRODIGAL.out.annotations )
            ch_versions       = ch_versions.mix(PYRODIGAL.out.versions)
            ch_annotation_faa = GUNZIP_PYRODIGAL_FAA.out.gunzip
            ch_annotation_fna = GUNZIP_PYRODIGAL_FNA.out.gunzip
            ch_annotation_gff = GUNZIP_PYRODIGAL_GFF.out.gunzip
            ch_annotation_gbk = Channel.empty() // Pyrodigal doesn't produce GBK

        }  else if ( params.annotation_tool == "prokka" ) {

            PROKKA ( fasta, [], [] )
            ch_versions       = ch_versions.mix(PROKKA.out.versions)
            ch_multiqc_files  = PROKKA.out.txt

            ch_annotation_faa = PROKKA.out.faa
            ch_annotation_fna = PROKKA.out.fna
            ch_annotation_gff = PROKKA.out.gff
            ch_annotation_gbk = PROKKA.out.gbk

        }   else if ( params.annotation_tool == "bakta" ) {

            // BAKTA prepare download
            if ( params.annotation_bakta_db_localpath ) {
                ch_bakta_db = Channel
                    .fromPath( params.annotation_bakta_db_localpath )
                    .first()
            } else {
                BAKTA_BAKTADBDOWNLOAD ( )
                ch_versions = ch_versions.mix( BAKTA_BAKTADBDOWNLOAD.out.versions )
                ch_bakta_db = ( BAKTA_BAKTADBDOWNLOAD.out.db )
            }

            BAKTA_BAKTA ( fasta, ch_bakta_db, [], [] )
            ch_versions       = ch_versions.mix(BAKTA_BAKTA.out.versions)
            ch_multiqc_files  = BAKTA_BAKTA.out.txt

            ch_annotation_faa = BAKTA_BAKTA.out.faa
            ch_annotation_fna = BAKTA_BAKTA.out.fna
            ch_annotation_gff = BAKTA_BAKTA.out.gff
            ch_annotation_gbk = BAKTA_BAKTA.out.gbff
        }

    emit:
    versions      = ch_versions
    multiqc_files = ch_multiqc_files
    faa           = ch_annotation_faa // [ [meta], path(faa) ]
    fna           = ch_annotation_fna // [ [meta], path(fna) ]
    gff           = ch_annotation_gff // [ [meta], path(gff) ]
    gbk           = ch_annotation_gbk // [ [meta], path(gbk) ]
}
