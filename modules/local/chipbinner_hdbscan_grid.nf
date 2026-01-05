process CHIPBINNER_HDBSCAN_GRID {
    tag "${group}"
    label 'process_high'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/02_chipbinner/${group}/clustering" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "conda-forge::python=3.11 conda-forge::numpy conda-forge::pandas conda-forge::scikit-learn conda-forge::hdbscan"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hdbscan:0.8.33--pyhdfd78af_0' :
        'biocontainers/hdbscan:0.8.33--pyhdfd78af_0' }"

    input:
    path matrix
    val samples_json
    val group
    val min_cluster_size
    val min_samples

    output:
    path "chipbinner.hdbscan_grid_summary.tsv", emit: grid_summary
    path "chipbinner.clusters.best.tsv"       , emit: clusters
    path "chipbinner.clusters.2clusters.tsv"  , emit: clusters_2
    path "chipbinner.clusters.3clusters.tsv"  , emit: clusters_3
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """\
    python - <<'PY'
    import json
    import csv

    samples = json.loads(r'''${samples_json}''')
    with open('chipbinner.samplesheet.csv', 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['sample_id', 'condition'])
        for row in samples:
            writer.writerow([row['sample_id'], row['condition']])
    PY

    python ${projectDir}/bin/hdbscan_grid.py \
        --matrix ${matrix} \
        --samplesheet chipbinner.samplesheet.csv \
        --min_cluster_size ${min_cluster_size} \
        --min_samples ${min_samples} \
        --summary chipbinner.hdbscan_grid_summary.tsv \
        --clusters_best chipbinner.clusters.best.tsv \
        --clusters_2 chipbinner.clusters.2clusters.tsv \
        --clusters_3 chipbinner.clusters.3clusters.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """.stripIndent()
}
