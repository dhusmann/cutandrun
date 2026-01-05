process CACHE_DIFFERENTIAL_GENE_BED {
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path gene_bed

    output:
    path "differential_manifest.gene_bed.bed", emit: bed
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cp ${gene_bed} differential_manifest.gene_bed.bed

    coreutils_version=\$(cat --version | head -n 1 | awk '{print $NF}')
    {
        echo "\"${task.process}\":"
        echo "    coreutils: $coreutils_version"
    } > versions.yml
    """
}
