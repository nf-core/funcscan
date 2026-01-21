process BAKTA_BAKTADBDOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bakta:1.11.4--pyhdfd78af_0'
        : 'biocontainers/bakta:1.11.4--pyhdfd78af_0'}"

    output:
    path "db*/", emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def db_type = args.tokenize()[args.tokenize().indexOf('--type') + 1]
    def url = db_type == 'light' ? "https://zenodo.org/records/14916843/files/db-light.tar.xz" : "https://zenodo.org/records/14916843/files/db.tar.xz"
    """
    wget ${url} 

    bakta_db \\
        install \\
        --db-file \$(find -name '*.xz')
        # ${args} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta_db --version) 2>&1 | cut -f '2' -d ' ')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo "bakta_db \\
        install \\
        --db-file *.xz
        # ${args} \\

    mkdir db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta_db --version) 2>&1 | cut -f '2' -d ' ')
    END_VERSIONS
    """
}
