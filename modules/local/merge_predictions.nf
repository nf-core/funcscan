process MERGE_PREDICTIONS {
    tag "$samplesheet"

    conda (params.enable_conda ? "conda-forge::python=3.8.6 conda-forge::pandas=1.0.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0' :
        'quay.io/biocontainers/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0' }"

    input:
    tuple val(meta), path(macrel), path(ampir), ...

    output:
    path output        , emit: predictions
    path 'versions.yml', emit: versions

    script:
    output = task.ext.prefix ?: 'predictions.tsv.gz'
    """
    merge_predictions.py \\
        --output $output \\
        --macrel $macrel \\
        --ampir $ampir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
