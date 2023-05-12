process DRAMP_DOWNLOAD {
    label 'process_single'

    conda "bioconda::ampcombi=0.1.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampcombi:0.1.7--pyhdfd78af_0':
        'biocontainers/ampcombi:0.1.7--pyhdfd78af_0' }"

    output:
    path "amp_ref_database/"    , emit: db
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/funcscan/bin/
    """
    mkdir amp_ref_database/
    ampcombi_download.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: \$(ampcombi --version | sed 's/ampcombi //')
    END_VERSIONS
    """
}
