process PEAK_QC_TABLE_REPORT {
    tag "$table_id"
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(table_id), val(header), val(rows)

    output:
    path "${table_id}.tsv", emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def lines = rows ? rows.join('\n') : ''
    """
    printf "%s\\n" "${header}" > ${table_id}.tsv
    if [ -n "${lines}" ]; then
        printf "%s\\n" "${lines}" >> ${table_id}.tsv
    fi

    coreutils_version=\$(cat --version | head -n 1 | awk '{print \$NF}')
    {
        echo "\\"${task.process}\\":"
        echo "    coreutils: \$coreutils_version"
    } > versions.yml
    """
}
