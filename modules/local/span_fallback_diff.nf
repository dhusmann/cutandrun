process SPAN_FALLBACK_DIFF {
    tag "${group}"
    label 'process_medium'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/03_span/${group}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "bioconda::bioconductor-diffbind bioconda::bioconductor-deseq2 bioconda::bioconductor-edger bioconda::bioconductor-genomicranges bioconda::bioconductor-genomicalignments bioconda::bioconductor-rsamtools conda-forge::r-base=4.2.3 conda-forge::r-optparse"
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
    val backend

    output:
    tuple val(group), path("span_fallback.differential.tsv")     , emit: results
    tuple val(group), path("span_fallback.differential.bed")     , emit: bed
    tuple val(group), path("span_fallback.significant.bed")      , emit: significant
    tuple val(group), path("span_fallback.significant_up.bed")   , emit: up
    tuple val(group), path("span_fallback.significant_down.bed") , emit: down
    tuple val(group), path("span_fallback.summary.tsv")          , emit: summary
    tuple val(group), path("span_fallback.readme.txt")           , emit: readme
    tuple val(group), path("span_fallback.union.bed")            , emit: union
    tuple val(group), path("span_fallback.counts.tsv")           , emit: counts
    path "00_manifests/span_fallback_union_manifest.tsv"         , emit: manifest
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript ${projectDir}/bin/span_fallback_de.R \\
        --samplesheet ${samplesheet} \\
        --treated ${treated} \\
        --control ${control} \\
        --group ${group} \\
        --fdr ${fdr} \\
        --backend ${backend} \\
        --use_spikein ${use_spikein} \\
        --prefix span_fallback

    mv span_fallback.results.tsv span_fallback.differential.tsv
    awk 'BEGIN{OFS="\\t"} NR>1 {print \$1,\$2,\$3}' span_fallback.differential.tsv > span_fallback.differential.bed

    mkdir -p 00_manifests
    peaks_list=\$(cat span_fallback.peaks.list 2>/dev/null || true)
    cat <<-END_MANIFEST > 00_manifests/span_fallback_union_manifest.tsv
    group\tunion_bed\tcount_matrix\tsource_peaks
    ${group}\t${params.outdir}/03_peak_calling/06_differential/03_span/${group}/span_fallback.union.bed\t${params.outdir}/03_peak_calling/06_differential/03_span/${group}/span_fallback.counts.tsv\t\${peaks_list}
    END_MANIFEST
    rm -f span_fallback.peaks.list

    cat <<-END_README > span_fallback.readme.txt
    SPAN fallback differential: SPAN peaks -> union peaks -> per-sample counts -> ${backend} differential.
    Log2FC is treated (${treated}) vs control (${control}). FDR threshold: ${fdr}.
    Spike-in scaling applied: ${use_spikein}.
    END_README

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deseq2: \$(Rscript -e 'packageVersion("DESeq2")' 2>/dev/null | tr -d '[]')
        r-base: \$(R --version 2>&1 | head -n 1 | sed -e 's/.*R version //; s/ .*//')
    END_VERSIONS
    """
}
