process CONTROL_POOLING_FALLBACKS_REPORT {
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    val records

    output:
    path "control_pooling_fallbacks.tsv", emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def header = "sample_id\tgroup\tcondition\tcontrol_group\tcontrol_condition\tcaller\treason"
    def lines = records ? records.collect { record ->
        [
            record.sample_id,
            record.group,
            record.condition,
            record.control_group,
            record.control_condition,
            record.caller,
            record.reason
        ].join('\t')
    }.join('\n') : ''
    """
    printf "%s\\n" "${header}" > control_pooling_fallbacks.tsv
    if [ -n "${lines}" ]; then
        printf "%s\\n" "${lines}" >> control_pooling_fallbacks.tsv
    fi

    coreutils_version=\$(cat --version | head -n 1 | awk '{print \$NF}')
    {
        echo "\\"${task.process}\\":"
        echo "    coreutils: \$coreutils_version"
    } > versions.yml
    """
}
