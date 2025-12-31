process NORMALISATION_FACTORS_REPORT {
    tag "$scope_id"
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(scope_id), val(records)

    output:
    path "${scope_id}.tsv", emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def header = "sample_id\tgroup\tcondition\treplicate\tspikein_reads\tscale_factor\tscope_id"
    def lines = records ? records.collect { record ->
        [
            record.sample_id,
            record.group,
            record.condition,
            record.replicate,
            record.spikein_reads,
            record.scale_factor,
            record.scope_id
        ].join('\t')
    }.join('\n') : ''
    """
    printf "%s\\n" "${header}" > ${scope_id}.tsv
    if [ -n "${lines}" ]; then
        printf "%s\\n" "${lines}" >> ${scope_id}.tsv
    fi

    coreutils_version=\$(cat --version | head -n 1 | awk '{print \$NF}')
    {
        echo "\\"${task.process}\\":"
        echo "    coreutils: \$coreutils_version"
    } > versions.yml
    """
}
