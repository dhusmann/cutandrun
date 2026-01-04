process CHIPBINNER_RUN {
    label 'process_chipbinner'

    conda "conda-forge::python=3.8.3"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    tuple val(group), path(records), path(chrom_sizes)
    val contrast
    val bin_size
    val windows_dir
    val use_input
    val pseudocount
    val grid_minpts
    val grid_minsamps
    val fdr
    val lfc
    val bootstrap
    val k_value
    val functional_db

    output:
    tuple val(group), path("chipbinner.samplesheet.csv"), emit: samplesheet
    tuple val(group), path("chipbinner.windows.bed"), emit: windows
    tuple val(group), path("chipbinner.bin_counts.tsv"), emit: counts
    tuple val(group), path("chipbinner.normalized_matrix.tsv"), emit: normalized
    tuple val(group), path("chipbinner.hdbscan_grid_summary.tsv"), emit: grid
    tuple val(group), path("chipbinner.clusters.tsv"), emit: clusters
    tuple val(group), path("chipbinner.differential.tsv"), emit: differential
    tuple val(group), path("chipbinner.summary.tsv"), emit: summary
    tuple val(group), path("plots"), emit: plots
    tuple val(group), path("enrichment"), emit: enrichment
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def windows_arg = windows_dir ? "--windows-dir ${windows_dir}" : ''
    def use_input_arg = use_input ? "--use-input" : ''
    def functional_arg = functional_db ? "--functional-db ${functional_db}" : ''
    """
    chipbinner_run.py \
        --samples ${records} \
        --group '${group}' \
        --contrast '${contrast}' \
        --chrom-sizes ${chrom_sizes} \
        --bin-size ${bin_size} \
        ${windows_arg} \
        ${use_input_arg} \
        --pseudocount ${pseudocount} \
        --grid-minpts '${grid_minpts}' \
        --grid-minsamps '${grid_minsamps}' \
        --fdr ${fdr} \
        --lfc ${lfc} \
        --bootstrap ${bootstrap} \
        --k-value ${k_value} \
        ${functional_arg} \
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o "([0-9]{1,}\\.)+[0-9]{1,}")
    END_VERSIONS
    """
}
