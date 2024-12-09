process INTERPROSCAN_DATABASE {
    tag "interproscan_database_download"
    //label 'process_medium'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/curl:7.80.0' :
        'biocontainers/curl:7.80.0' }"

    input:
    val (database_url)

    output:
    path("interproscan_db/*")    , emit: db
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir interproscan_db/

    curl -L ${database_url} -o interproscan_db/interproscan-5.59-91.0-64-bit.tar.gz
    tar -xzf interproscan_db/interproscan-5.59-91.0-64-bit.tar.gz -C interproscan_db/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version 2>&1 | sed -n '1s/tar (busybox) //p')
        curl: "\$(curl --version 2>&1 | sed -n '1s/^curl \\([0-9.]*\\).*/\\1/p')"
    END_VERSIONS
    """
}
