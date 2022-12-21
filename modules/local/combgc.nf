process COMBGC {
    tag "comBGC"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75' :
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    path input_antismash
    path input_deepbgc
    path input_gecco

    output:
    path 'results/combgc/combgc_summary.tsv' , emit: tsv

    script: // This script is bundled with the pipeline, in nf-core/funcscan/bin/
    def arg_antismash = input_antismash ? "--amp_database $input_antismash": ""
    def arg_deepbgc = input_deepbgc ? "--amp_database $input_deepbgc": ""
    def arg_gecco = input_gecco ? "--amp_database $input_gecco": ""
    """
    comBGC.py \\
        $arg_antismash \\
        $arg_deepbgc \\
        $arg_gecco
    """
}
