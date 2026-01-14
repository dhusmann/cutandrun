process PUBLISH_DIFFERENTIAL_MANIFESTS {
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path samples_manifest
    path peaks_manifest

    output:
    path "differential_manifest.samples.tsv", emit: samples
    path "differential_manifest.peaks.tsv", emit: peaks
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samples_src=\$(readlink -f ${samples_manifest})
    samples_dst=\$(readlink -f differential_manifest.samples.tsv)
    if [[ "\$samples_src" != "\$samples_dst" ]]; then
        cp ${samples_manifest} differential_manifest.samples.tsv
    fi

    peaks_src=\$(readlink -f ${peaks_manifest})
    peaks_dst=\$(readlink -f differential_manifest.peaks.tsv)
    if [[ "\$peaks_src" != "\$peaks_dst" ]]; then
        cp ${peaks_manifest} differential_manifest.peaks.tsv
    fi

    coreutils_version=\$(cat --version | head -n 1 | awk '{print \$NF}')
    {
        echo "\"${task.process}\":"
        echo "    coreutils: \$coreutils_version"
    } > versions.yml
    """
}
