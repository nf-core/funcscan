process COMBGC {
    tag "comBGC"

    conda "conda-forge::python=3.9.0 conda-forge::biopython=1.79 conda-forge::pandas=1.3.5" // Using versions from available mulled containers. Most recent would instead be: conda-forge::python=3.11.0 conda-forge::biopython=1.80 conda-forge::pandas=1.5.2
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' :
        'quay.io/biocontainers/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' }"

    input:
    path antismash_dir
    path deepbgc_dir
    path gecco_dir
    path out_dir

    output:
    path "${out_dir}/combgc_summary.tsv" , emit: tsv

    script: // This script is bundled with the pipeline, in nf-core/funcscan/bin/
    def antismash_dir = antismash_dir ? "--antismash $antismash_dir": ""
    def deepbgc_dir = deepbgc_dir ? "--deepbgc $deepbgc_dir": ""
    def gecco_dir = gecco_dir ? "--gecco $gecco_dir": ""
    def out_dir = out_dir ? "--outdir $out_dir": ""
    """
    comBGC.py \\
        $antismash_dir \\
        $deepbgc_dir \\
        $gecco_dir \\
        $out_dir
    """
}
