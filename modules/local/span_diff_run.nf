process SPAN_DIFF_RUN {
    label 'process_span'

    conda "conda-forge::python=3.8.3" 
    container "quay.io/biocontainers/python:3.8.3"

    input:
    tuple val(group), path(records), path(peaks)
    val contrast
    val mode
    val fdr
    val gap
    val bin_size
    val fallback_backend
    path omnipeaks_jar
    val java_heap

    output:
    tuple val(group), path("span.differential.tsv"), emit: differential
    tuple val(group), path("span.differential.peaks.bed"), emit: peaks_bed
    tuple val(group), path("span.up.bed"), emit: up
    tuple val(group), path("span.down.bed"), emit: down
    tuple val(group), path("span.summary.tsv"), emit: summary
    tuple val(group), path("span.mode.txt"), emit: mode_out
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    span_diff.py \
        --jar ${omnipeaks_jar} \
        --mode ${mode} \
        --contrast '${contrast}' \
        --group '${group}' \
        --samples ${records} \
        --peaks ${peaks} \
        --bin ${bin_size} \
        --gap ${gap} \
        --fdr ${fdr} \
        --fallback-backend ${fallback_backend} \
        --java-heap ${java_heap} \
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o "([0-9]{1,}\\.)+[0-9]{1,}")
    END_VERSIONS
    """
}
