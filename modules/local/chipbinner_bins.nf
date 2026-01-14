process CHIPBINNER_BINS {
    tag "${group}.${bin_size}"
    label 'process_high'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/02_chipbinner/${group}/bins" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    path chrom_sizes
    val bin_size
    val group
    val windows_dir
    val blacklist

    output:
    tuple val(group), path("*.windows.bed"), emit: bins
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${group}.${bin_size}"
    def windows_arg = windows_dir ? "${windows_dir}" : ''
    def blacklist_arg = blacklist ? "${blacklist}" : ''
    """
    set -euo pipefail
    if [ -n "${windows_arg}" ]; then
        windows_file=\$(ls ${windows_arg}/*${bin_size}*.bed 2>/dev/null | head -n 1 || true)
        if [ -z "\$windows_file" ]; then
            echo "No window file found in ${windows_arg} for bin size ${bin_size}" >&2
            exit 1
        fi
        cp "\$windows_file" ${prefix}.windows.bed
    else
        bedtools makewindows -g ${chrom_sizes} -w ${bin_size} > ${prefix}.windows.bed
    fi

    if [ -n "${blacklist_arg}" ]; then
        bedtools subtract -a ${prefix}.windows.bed -b ${blacklist_arg} > ${prefix}.windows.filtered.bed
        mv ${prefix}.windows.filtered.bed ${prefix}.windows.bed
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
