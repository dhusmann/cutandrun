process DIFFERENTIAL_SUMMARY {
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path design_manifest
    path summary_files
    path summary_header
    path design_header

    output:
    path "differential_summary_mqc.tsv", emit: summary
    path "differential_design_mqc.tsv", emit: design
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    differential_summary.py \
        --design ${design_manifest} \
        --summaries ${summary_files} \
        --out differential_summary.tsv

    cat ${summary_header} differential_summary.tsv > differential_summary_mqc.tsv
    cat ${design_header} ${design_manifest} > differential_design_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o "([0-9]{1,}\\.)+[0-9]{1,}")
    END_VERSIONS
    """
}
