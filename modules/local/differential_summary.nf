process DIFFERENTIAL_SUMMARY {
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path design_manifest
    path summary_files
    val manifest_only
    path summary_header
    path design_header
    path chipbinner_header
    path span_header

    output:
    path "differential_summary_mqc.tsv", emit: summary
    path "differential_design_mqc.tsv", emit: design
    path "chipbinner_summary_mqc.tsv", emit: chipbinner_summary
    path "span_summary_mqc.tsv", emit: span_summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def summaries_arg = summary_files instanceof List ? summary_files.join(' ') : summary_files
    def manifest_only_flag = manifest_only ? '--manifest-only' : ''
    """
    differential_summary.py \
        --design ${design_manifest} \
        --summaries ${summaries_arg} \
        ${manifest_only_flag} \
        --out differential_summary.tsv

    cat ${summary_header} differential_summary.tsv > differential_summary_mqc.tsv
    cat ${design_header} ${design_manifest} > differential_design_mqc.tsv
    awk 'BEGIN{FS=OFS="\\t"} NR==1{print; next} \$1=="chipbinner"' differential_summary.tsv > chipbinner_summary.tsv
    cat ${chipbinner_header} chipbinner_summary.tsv > chipbinner_summary_mqc.tsv
    awk 'BEGIN{FS=OFS="\\t"} NR==1{print; next} \$1=="span"' differential_summary.tsv > span_summary.tsv
    cat ${span_header} span_summary.tsv > span_summary_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o "([0-9]{1,}\\.)+[0-9]{1,}")
    END_VERSIONS
    """
}
