process RECORDS_TO_TSV {
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    val key1
    val key2
    val records_json
    val header
    val filename

    output:
    tuple val(key1), val(key2), path("${filename}"), emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    records_to_tsv.py \
        --records '${records_json}' \
        --header '${header}' \
        --out '${filename}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o "([0-9]{1,}\\.)+[0-9]{1,}")
    END_VERSIONS
    """
}
