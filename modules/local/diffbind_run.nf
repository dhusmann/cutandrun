process DIFFBIND_RUN {
    label 'process_diffbind'

    conda "conda-forge::r-base=4.2.3 bioconda::bioconductor-diffbind conda-forge::r-jsonlite conda-forge::r-yaml"
    container "quay.io/biocontainers/bioconductor-diffbind:3.14.0--r42_0"

    input:
    tuple val(group), val(caller), path(records)
    val contrast
    val use_spikein
    val fdr
    val lfc
    val min_overlap
    val backend
    val recenter
    val summits
    val norm_method
    val extra_params
    val export_sheets

    output:
    tuple val(group), val(caller), path("diffbind.results.tsv"), emit: results
    tuple val(group), val(caller), path("diffbind.significant.bed"), emit: bed
    tuple val(group), val(caller), path("diffbind.significant_up.bed"), emit: bed_up
    tuple val(group), val(caller), path("diffbind.significant_down.bed"), emit: bed_down
    tuple val(group), val(caller), path("diffbind.summary.tsv"), emit: summary
    tuple val(group), val(caller), path("diffbind.samplesheet.csv"), emit: samplesheet
    tuple val(group), val(caller), path("diffbind.normalization_factors.tsv"), emit: norm_factors_out
    tuple val(group), val(caller), path("diffbind.dba.rds"), emit: dba
    tuple val(group), val(caller), path("plots"), emit: plots
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def extra_arg = extra_params ? "--extra_params ${extra_params}" : ''
    """
    diffbind_run.R \
        --records ${records} \
        --outdir . \
        --contrast '${contrast}' \
        --caller '${caller}' \
        --group '${group}' \
        --use_spikein ${use_spikein} \
        --fdr ${fdr} \
        --lfc ${lfc} \
        --min_overlap ${min_overlap} \
        --backend ${backend} \
        --recenter ${recenter} \
        --summits ${summits} \
        --norm_method ${norm_method} \
        ${extra_arg} \
        --export_sheets ${export_sheets}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1 | sed 's/.* //')
    END_VERSIONS
    """
}
