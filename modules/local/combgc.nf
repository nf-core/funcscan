process COMBGC {
    tag "comBGC"

    conda (params.enable_conda ? "conda-forge::python=3.10.2 conda-forge::biopython=1.79 conda-forge::pandas=1.5.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75' :
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    path antismash_dir
    path deepbgc_dir
    path gecco_dir
    path out_dir

    output:
    path "${out_dir}/combgc_summary.tsv" , emit: tsv

    script: // This script is bundled with the pipeline, in nf-core/funcscan/bin/
    def arg_antismash = antismash_dir ? "--antismash $antismash_dir": ""
    def arg_deepbgc = deepbgc_dir ? "--deepbgc $deepbgc_dir": ""
    def arg_gecco = gecco_dir ? "--gecco $gecco_dir": ""
    def arg_outdir = out_dir ? "--outdir $out_dir": ""
    """
    comBGC.py \\
        $arg_antismash \\
        $arg_deepbgc \\
        $arg_gecco \\
        $arg_outdir
    """
}
