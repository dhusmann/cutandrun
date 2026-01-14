process DIFFERENTIAL_MANIFESTS {
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path samples_raw
    path peaks_raw
    path ms_coeffs
    val normalisation_mode

    output:
    path "differential_manifest.samples.tsv", emit: samples
    path "differential_manifest.peaks.tsv", emit: peaks
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ms_arg = ms_coeffs ? "--ms-coeffs ${ms_coeffs}" : ''
    """
    differential_manifests.py \
        --samples ${samples_raw} \
        --peaks ${peaks_raw} \
        --normalisation-mode ${normalisation_mode} \
        ${ms_arg} \
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o "([0-9]{1,}\\.)+[0-9]{1,}")
    END_VERSIONS
    """
}
