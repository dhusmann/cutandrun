process SPAN_FALLBACK_DIFF {
    tag "${group}"
    label 'process_medium'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/03_span/${group}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "bioconda::bioconductor-diffbind bioconda::bioconductor-deseq2 bioconda::bioconductor-edger bioconda::bioconductor-genomicranges conda-forge::r-base=4.2.3 conda-forge::r-optparse conda-forge::r-ggplot2 conda-forge::r-jsonlite conda-forge::r-yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-diffbind:3.14.0--r42hdfd78af_0' :
        'biocontainers/bioconductor-diffbind:3.14.0--r42hdfd78af_0' }"

    input:
    path samplesheet
    val treated
    val control
    val group
    val use_spikein
    val fdr
    val lfc
    val min_overlap
    val summits
    val backend

    output:
    tuple val(group), path("span_fallback.differential.tsv")     , emit: results
    tuple val(group), path("span_fallback.differential.bed")     , emit: bed
    tuple val(group), path("span_fallback.significant.bed")      , emit: significant
    tuple val(group), path("span_fallback.significant_up.bed")   , emit: up
    tuple val(group), path("span_fallback.significant_down.bed") , emit: down
    tuple val(group), path("span_fallback.summary.tsv")          , emit: summary
    tuple val(group), path("span_fallback.readme.txt")           , emit: readme
    tuple val(group), path("plots")                              , emit: plots
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p plots
    Rscript ${projectDir}/bin/diffbind_run.R \
        --samplesheet ${samplesheet} \
        --treated ${treated} \
        --control ${control} \
        --group ${group} \
        --caller span_fallback \
        --fdr ${fdr} \
        --lfc ${lfc} \
        --min_overlap ${min_overlap} \
        --summits ${summits} \
        --backend ${backend} \
        --use_spikein ${use_spikein} \
        --prefix span_fallback

    mv span_fallback.results.tsv span_fallback.differential.tsv
    cp span_fallback.significant.bed span_fallback.differential.bed

    cat <<-END_README > span_fallback.readme.txt
    SPAN fallback differential: SPAN peaks + DiffBind/DESeq2/edgeR analysis.
    This fallback runs DiffBind on SPAN-derived peak sets to approximate differential enrichment
    when native OmniPeaks differential mode is unavailable.
    END_README

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diffbind: \$(Rscript -e 'packageVersion("DiffBind")' 2>/dev/null | tr -d '[]')
        r-base: \$(R --version 2>&1 | head -n 1 | sed -e 's/.*R version //; s/ .*//')
    END_VERSIONS
    """
}
