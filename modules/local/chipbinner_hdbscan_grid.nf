process CHIPBINNER_HDBSCAN_GRID {
    tag "${group}"
    label 'process_high'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/02_chipbinner/${group}/clustering" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "conda-forge::python=3.11 conda-forge::numpy conda-forge::pandas conda-forge::scikit-learn conda-forge::hdbscan"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path matrix
    val group
    val min_cluster_size
    val min_samples

    output:
    path "chipbinner.hdbscan_grid_summary.tsv", emit: grid_summary
    path "chipbinner.clusters.best.tsv"       , emit: clusters
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """\
    python - <<'PY'
    import importlib.util
    import subprocess
    import sys

    required = {
        "numpy": "numpy",
        "pandas": "pandas",
        "sklearn": "scikit-learn",
        "hdbscan": "hdbscan",
    }
    missing = [pkg for mod, pkg in required.items() if importlib.util.find_spec(mod) is None]
    if missing:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "--no-cache-dir"] + missing)
    PY

    python ${projectDir}/bin/hdbscan_grid.py \
        --matrix ${matrix} \
        --min_cluster_size ${min_cluster_size} \
        --min_samples ${min_samples} \
        --summary chipbinner.hdbscan_grid_summary.tsv \
        --clusters chipbinner.clusters.best.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """.stripIndent()
}
