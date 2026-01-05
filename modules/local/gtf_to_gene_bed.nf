process GTF_TO_GENE_BED {
    tag "$gtf"
    label 'process_low'

    conda "conda-forge::python=3.11.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.6' :
        'quay.io/biocontainers/python:3.11.6' }"

    input:
    path gtf

    output:
    path "gene_annotation.bed", emit: bed
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gtf_to_gene_bed.py \\
        --gtf $gtf \\
        --out gene_annotation.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print $NF}')
    END_VERSIONS
    """
}
