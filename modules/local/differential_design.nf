process DIFFERENTIAL_DESIGN {
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path samples_manifest
    path peaks_manifest
    val contrast
    val min_reps
    val allow_partial
    val run_diffbind
    val run_chipbinner
    val run_span
    val groups_allow
    val callers_allow

    output:
    path "differential_manifest.design.tsv", emit: design
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def allow_flag = allow_partial ? '--allow-partial' : ''
    def diffbind_flag = run_diffbind ? '--run-diffbind' : ''
    def chipbinner_flag = run_chipbinner ? '--run-chipbinner' : ''
    def span_flag = run_span ? '--run-span' : ''
    def groups_arg = groups_allow ? "--groups ${groups_allow}" : ''
    def callers_arg = callers_allow ? "--callers ${callers_allow}" : ''
    """
    differential_design.py \
        --samples ${samples_manifest} \
        --peaks ${peaks_manifest} \
        --contrast '${contrast}' \
        --min-replicates ${min_reps} \
        ${allow_flag} \
        ${diffbind_flag} \
        ${chipbinner_flag} \
        ${span_flag} \
        ${groups_arg} \
        ${callers_arg} \
        --out differential_manifest.design.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o "([0-9]{1,}\\.)+[0-9]{1,}")
    END_VERSIONS
    """
}
