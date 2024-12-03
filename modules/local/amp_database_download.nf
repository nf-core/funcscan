process AMP_DATABASE_DOWNLOAD {
    label 'process_single'

    conda "bioconda::ampcombi=2.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampcombi:2.0.1--pyhdfd78af_0':
        'biocontainers/ampcombi:2.0.1--pyhdfd78af_0' }"

    input:
    val database_id

    output:
    path "amp_${database_id}_database"  , emit: db
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/funcscan/bin/
    """
    ampcombi_download.py \\
        --database_id $database_id \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: \$(ampcombi --version | sed 's/ampcombi //')
    END_VERSIONS
    """
}
