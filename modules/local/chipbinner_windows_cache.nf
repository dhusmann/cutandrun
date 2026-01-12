process CHIPBINNER_WINDOWS_CACHE {
    label 'process_single'

    conda "conda-forge::python=3.8.3 bioconda::bedtools=2.31.1"
    container workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        null :
        'quay.io/biocontainers/biocontainers:1.2.0--py38_0'

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

    stub:
    """
    genome_id=\$(basename "${chrom_sizes}")
    if [[ "\$genome_id" == *.sizes ]]; then
        genome_id="\${genome_id%.sizes}"
    else
        genome_id="\${genome_id%.*}"
    fi
    if [[ -z "\$genome_id" ]]; then
        genome_id="genome"
    fi

    blacklist_hash="none"
    if [[ -n "${blacklist}" && -f "${blacklist}" ]]; then
        if command -v sha256sum >/dev/null 2>&1; then
            blacklist_hash=\$(sha256sum "${blacklist}" | awk '{print substr(\$1,1,12)}')
        fi
    fi

    out_name="windows.\${genome_id}.${bin_size}.\${blacklist_hash}.bed"
    printf "chrStub\\t0\\t${bin_size}\\n" > "\${out_name}"

    {
        printf "genome_id\\tbin_size\\tblacklist_hash\\twindows_source\\twindows_path\\n"
        printf "%s\\t%s\\t%s\\t%s\\t%s\\n" "\$genome_id" "${bin_size}" "\$blacklist_hash" "${windows_dir ?: 'generated'}" "\${out_name}"
    } > chipbinner_windows_meta.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
        bedtools: "stub"
    END_VERSIONS
    """
}
