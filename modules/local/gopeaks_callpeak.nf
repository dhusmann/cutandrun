process GOPEAKS_CALLPEAK {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gopeaks=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gopeaks:1.0.0--h047eeb3_3' :
        'biocontainers/gopeaks:1.0.0--h047eeb3_3' }"

    input:
    tuple val(meta), path(bam)
    val   broad

    output:
    tuple val(meta), path("*.bed"), emit: peaks
    tuple val(meta), path("*_gopeaks.json"), emit: json
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.caller}"
    def broad_flag = broad ? '--broad' : ''
    """
    gopeaks -b $bam -o ${prefix} ${broad_flag} ${args}

    if [[ -f ${prefix}_peaks.bed ]]; then
        mv ${prefix}_peaks.bed ${prefix}.bed
    fi

    if [[ -f ${prefix}.json ]]; then
        mv ${prefix}.json ${prefix}_gopeaks.json
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gopeaks: \$(gopeaks --version 2>&1 | head -n 1 | sed -e 's/.*gopeaks //g' || echo "unknown")
    END_VERSIONS
    """
}
