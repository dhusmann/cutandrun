process SPAN_POOLING_MANIFEST {
    tag "${group}"
    label 'process_single'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/03_span/${group}/00_manifests" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    val group
    val records

    output:
    tuple val(group), path("span_pooling.tsv"), emit: manifest
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def header = "group\tcondition\tpooled_bam\tpooled_bai\tinput_bams\tpooling_reason\tcreated_by"
    def lines = records ? records.collect { record ->
        [
            record.group,
            record.condition,
            record.pooled_bam ?: 'NA',
            record.pooled_bai ?: 'NA',
            record.input_bams ?: '',
            record.pooling_reason ?: '',
            record.created_by ?: 'native'
        ].join('\t')
    }.join('\n') : ''
    """
    printf "%s\\n" "${header}" > span_pooling.tsv
    if [ -n "${lines}" ]; then
        printf "%s\\n" "${lines}" >> span_pooling.tsv
    fi

    coreutils_version=\$(cat --version | head -n 1 | awk '{print \$NF}')
    {
        echo "\\"${task.process}\\":"
        echo "    coreutils: \$coreutils_version"
    } > versions.yml
    """
}
