process PEAK_COUNTS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(bed)
    path  peak_counts_header

    output:
    tuple val(meta), path("*_peak_count.txt"), emit: count_value
    tuple val(meta), path("*mqc.tsv"), emit: count_mqc
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    PEAK_COUNT=\$(cat ${bed} | wc -l | awk '{print \$1}')
    printf "%s\\n" "\$PEAK_COUNT" > ${prefix}_peak_count.txt
    printf "Peak Count\\t%s\\n" "\$PEAK_COUNT" | cat $peak_counts_header - > ${prefix}_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
