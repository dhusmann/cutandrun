process CHIPBINNER_WINDOWS_CACHE {
    label 'process_single'

    conda "conda-forge::python=3.8.3 bioconda::bedtools=2.31.1"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path chrom_sizes
    val bin_size
    val windows_dir
    val blacklist

    output:
    path "windows.*.bed", emit: windows
    path "chipbinner_windows_meta.tsv", emit: meta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def windows_dir_arg = windows_dir ? "--windows-dir ${windows_dir}" : ''
    def blacklist_arg = blacklist ? "--blacklist ${blacklist}" : ''
    """
    chipbinner_windows_cache.py \
        --chrom-sizes ${chrom_sizes} \
        --bin-size ${bin_size} \
        ${windows_dir_arg} \
        ${blacklist_arg} \
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o "([0-9]{1,}\\.)+[0-9]{1,}")
        bedtools: \$(bedtools --version | awk '{print \$2}')
    END_VERSIONS
    """
}
