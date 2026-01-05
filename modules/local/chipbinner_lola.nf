process CHIPBINNER_LOLA {
    tag "${group}.${label}"
    label 'process_medium'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/02_chipbinner/${group}/enrichment" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "bioconda::bioconductor-lola=1.28.0 conda-forge::r-ggplot2 conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-lola:1.28.0--r42hdfd78af_0' :
        'biocontainers/bioconductor-lola:1.28.0--r42hdfd78af_0' }"

    input:
    tuple val(group), val(label), path(bed)
    path universe
    path lola_db

    output:
    tuple val(group), val(label), path("${label}.lola.tsv"), emit: tsv
    tuple val(group), val(label), path("${label}.lola_plot.pdf"), optional: true, emit: plot
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript ${projectDir}/bin/chipbinner_lola.R \
        --bed ${bed} \
        --universe ${universe} \
        --db ${lola_db} \
        --label ${label} \
        --out_tsv ${label}.lola.tsv \
        --out_pdf ${label}.lola_plot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | head -n 1 | sed -e 's/.*R version //; s/ .*//')
    END_VERSIONS
    """
}
