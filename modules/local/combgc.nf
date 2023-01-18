process COMBGC {
    tag "$meta.id"

    conda "conda-forge::python=3.9.0 conda-forge::biopython=1.79 conda-forge::pandas=1.3.5" // Using versions from available mulled containers. Most recent would instead be: conda-forge::python=3.11.0 conda-forge::biopython=1.80 conda-forge::pandas=1.5.2
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' :
        'quay.io/biocontainers/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' }"

    input:
    tuple val(meta), path(antismash_dir, stageAs:"antismash/*"), path(deepbgc_dir, stageAs:"deepbgc/*"), path(gecco_dir, stageAs:"gecco/*")

    output:
    path "${prefix}/combgc_summary.tsv" , emit: tsv

    script: // This script is bundled with the pipeline, in nf-core/funcscan/bin/
    prefix = task.ext.prefix ?: "${meta.id}"
    prep_antismash = antismash_dir ? "mkdir -p antismash_in/${prefix}; mv antismash/* antismash_in/${prefix}/" : ""
    prep_deepbgc = deepbgc_dir ? "mkdir -p deepbgc_in/${prefix}; mv deepbgc/* deepbgc_in/${prefix}/" : ""
    prep_gecco = gecco_dir ? "mkdir -p gecco_in/${prefix}; mv gecco/* gecco_in/${prefix}/" : ""

    antismash_arg = antismash_dir ? "-a antismash_in/${prefix}" : ""
    deepbgc_arg = deepbgc_dir ? "-d deepbgc_in/${prefix}" : ""
    gecco_arg = gecco_dir ? "-g gecco_in/${prefix}" : ""
    """
    $prep_antismash
    $prep_deepbgc
    $prep_gecco

    comBGC.py \\
        $antismash_arg \\
        $deepbgc_arg \\
        $gecco_arg \\
        -o $prefix
    """
}
