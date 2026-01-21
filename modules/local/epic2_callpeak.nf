process EPIC2_CALLPEAK {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::epic2=0.0.54"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epic2:0.0.54--py312hdcc493e_0' :
        'biocontainers/epic2:0.0.54--py312hdcc493e_0' }"

    input:
    tuple val(meta), path(treatment_bam), path(control_bam)
    val   epic2_genome
    val   window_size
    val   gaps_allowed
    val   fdr

    output:
    tuple val(meta), path("*.peaks"), emit: peaks
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.caller}"
    """
    epic2 \
        --treatment ${treatment_bam} \
        --control ${control_bam} \
        --genome ${epic2_genome} \
        --bin-size ${window_size} \
        --gaps-allowed ${gaps_allowed} \
        --false-discovery-rate-cutoff ${fdr} \
        --output ${prefix}.peaks \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        epic2: \$(epic2 --version 2>&1 | grep -Eo '[0-9]+(\\.[0-9]+)+' | head -n 1 || echo "unknown")
    END_VERSIONS
    """
}
