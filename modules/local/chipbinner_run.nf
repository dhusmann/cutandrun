process CHIPBINNER_RUN {
    label 'process_chipbinner'

    conda "conda-forge::python=3.8.3 conda-forge::numpy=1.24.4 conda-forge::pandas=2.0.3 conda-forge::scikit-learn=1.3.2 conda-forge::matplotlib=3.7.5 conda-forge::seaborn=0.12.2 conda-forge::hdbscan=0.8.33 conda-forge::scipy=1.10.1 bioconda::bedtools=2.31.1 bioconda::samtools=1.17 conda-forge::r-base=4.2.3 bioconda::bioconductor-rots=1.20.0"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    tuple val(group), path(records), path(chrom_sizes), path(windows)
    val contrast
    val bin_size
    val windows_dir
    val blacklist
    val use_input
    val use_spikein
    val pseudocount
    val grid_minpts
    val grid_minsamps
    val fdr
    val lfc
    val bootstrap
    val k_value
    val functional_db
    val allow_partial

    output:
    tuple val(group), path("chipbinner.samplesheet.csv"), emit: samplesheet, optional: true
    tuple val(group), path("chipbinner.windows.bed"), emit: windows, optional: true
    tuple val(group), path("chipbinner.bin_counts.tsv"), emit: counts, optional: true
    tuple val(group), path("chipbinner.normalized_matrix.tsv"), emit: normalized, optional: true
    tuple val(group), path("chipbinner.normalization_factors.tsv"), emit: norm_factors, optional: true
    tuple val(group), path("chipbinner.hdbscan_grid_summary.tsv"), emit: grid, optional: true
    tuple val(group), path("hdbscan_grid"), emit: grid_outputs, optional: true
    tuple val(group), path("chipbinner.clusters.tsv"), emit: clusters, optional: true
    tuple val(group), path("chipbinner.clusters.best.tsv"), emit: clusters_best, optional: true
    tuple val(group), path("chipbinner.clusters.2cluster.tsv"), emit: clusters_two, optional: true
    tuple val(group), path("chipbinner.clusters.3cluster.tsv"), emit: clusters_three, optional: true
    tuple val(group), path("chipbinner.differential.tsv"), emit: differential, optional: true
    tuple val(group), path("chipbinner.summary.tsv"), emit: summary
    tuple val(group), path("chipbinner.error.txt"), emit: error, optional: true
    tuple val(group), path("plots"), emit: plots, optional: true
    tuple val(group), path("enrichment"), emit: enrichment, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def windows_arg = windows ? "--windows ${windows}" : ''
    def windows_dir_arg = (!windows && windows_dir) ? "--windows-dir ${windows_dir}" : ''
    def blacklist_arg = (!windows && blacklist) ? "--blacklist ${blacklist}" : ''
    def use_input_arg = use_input ? "--use-input" : ''
    def use_spikein_arg = use_spikein ? "--use-spikein" : ''
    def functional_arg = functional_db ? "--functional-db ${functional_db}" : ''
    def allow_partial_arg = allow_partial ? "--allow-partial" : ''
    """
    chipbinner_run.py \
        --samples ${records} \
        --group '${group}' \
        --contrast '${contrast}' \
        --chrom-sizes ${chrom_sizes} \
        --bin-size ${bin_size} \
        ${windows_arg} \
        ${windows_dir_arg} \
        ${blacklist_arg} \
        ${use_input_arg} \
        ${use_spikein_arg} \
        --pseudocount ${pseudocount} \
        --grid-minpts '${grid_minpts}' \
        --grid-minsamps '${grid_minsamps}' \
        --fdr ${fdr} \
        --lfc ${lfc} \
        --bootstrap ${bootstrap} \
        --k-value ${k_value} \
        ${functional_arg} \
        ${allow_partial_arg} \
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o "([0-9]{1,}\\.)+[0-9]{1,}")
    END_VERSIONS
    """
}
