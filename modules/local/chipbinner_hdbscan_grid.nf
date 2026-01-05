process CHIPBINNER_HDBSCAN_GRID {
    tag "${group}"
    label 'process_high'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/02_chipbinner/${group}/clustering" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "conda-forge::python=3.11 conda-forge::numpy conda-forge::hdbscan"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hdbscan:0.8.33--pyhdfd78af_0' :
        'biocontainers/hdbscan:0.8.33--pyhdfd78af_0' }"

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
