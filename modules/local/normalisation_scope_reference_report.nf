process NORMALISATION_SCOPE_REFERENCE_REPORT {
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    val records

    output:
    path "normalisation_scope_reference.tsv", emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def header = "scope_id\treference_reads"
    def lines = records ? records.collect { record ->
        [
            record.scope_id,
            record.reference_reads
        ].join('\t')
    }.join('\n') : ''
    """
    printf "%s\\n" "${header}" > normalisation_scope_reference.tsv
    if [ -n "${lines}" ]; then
        printf "%s\\n" "${lines}" >> normalisation_scope_reference.tsv
    fi

    coreutils_version=\$(cat --version | head -n 1 | awk '{print \$NF}')
    {
        echo "\\"${task.process}\\":"
        echo "    coreutils: \$coreutils_version"
    } > versions.yml
    """
}
