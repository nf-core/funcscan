process MERGE_TAXONOMY_COMBGC {
    label 'process_medium'

    conda "conda-forge::python=3.11.0 conda-forge::biopython=1.80 conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' :
        'biocontainers/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' }"

    input:
    path(combgc_df)
    path(taxa_list)

    output:
    path "combgc_complete_summary_taxonomy.tsv" , emit: tsv
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/funcscan/bin/
    """
    merge_taxonomy.py \\
        combgc_taxa \\
        --combgc $combgc_df \\
        --taxonomy $taxa_list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge_taxonomy: \$(merge_taxonomy.py --version | sed 's/merge_taxonomy //g')
    END_VERSIONS
    """
}
