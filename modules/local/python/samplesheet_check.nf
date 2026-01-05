process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path samplesheet

    output:
    path '*.csv'        , emit: csv
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def allow_cross = params.allow_cross_condition_controls ? '--allow-cross-condition-controls' : ''
    """
    check_samplesheet.py $samplesheet samplesheet.valid.csv $params.use_control $allow_cross

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o \"([0-9]{1,}\\.)+[0-9]{1,}\")
    END_VERSIONS
    """
}
