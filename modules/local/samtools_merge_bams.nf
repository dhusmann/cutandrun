process SAMTOOLS_MERGE_BAMS {
    tag "${meta.group}.${meta.role}"
    label 'process_medium'

    conda "bioconda::samtools=1.19.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("${meta.group}.${meta.role}.merged.bam"), path("${meta.group}.${meta.role}.merged.bam.bai"), emit: merged
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.group}.${meta.role}.merged"
    """
    samtools merge -f ${prefix}.bam ${bams}
    samtools index ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -n 1 | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
