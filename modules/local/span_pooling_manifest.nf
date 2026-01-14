process SPAN_POOLING_MANIFEST {
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path pooling_files
    val do_run

    output:
    path "span_diff_target_pooling.tsv", emit: manifest
    path "versions.yml", emit: versions

    when:
    (task.ext.when == null || task.ext.when) && do_run

    script:
    def inputs_arg = pooling_files instanceof List ? pooling_files.join(' ') : pooling_files
    """
    span_pooling_manifest.py \
        --inputs ${inputs_arg} \
        --out span_diff_target_pooling.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o "([0-9]{1,}\\\\.)+[0-9]{1,}")
    END_VERSIONS
    """
}
