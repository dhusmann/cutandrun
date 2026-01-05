process SPAN_DIFF_RUN {
    label { params.span_diff_mode == 'fallback' ? 'SPAN_FALLBACK' : 'SPAN_NATIVE' }

    conda "conda-forge::python=3.8.3 conda-forge::openjdk=21.0.2 bioconda::samtools=1.16.1 bioconda::bedtools=2.31.0 conda-forge::r-base=4.2.3 bioconda::bioconductor-deseq2 bioconda::bioconductor-edger"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    tuple val(group), path(samples_manifest), path(peaks_manifest), path(chrom_sizes)
    val contrast
    val mode
    val fdr
    val gap
    val bin_size
    val fallback_backend
    val use_spikein
    val allow_partial
    val caller_priority
    path omnipeaks_jar
    val java_heap

    output:
    tuple val(group), path("span.differential.tsv"), emit: differential
    tuple val(group), path("span.differential.peaks.bed"), emit: peaks_bed
    tuple val(group), path("span.up.bed"), emit: up
    tuple val(group), path("span.down.bed"), emit: down
    tuple val(group), path("span.summary.tsv"), emit: summary
    tuple val(group), path("span.mode.txt"), emit: mode_out
    tuple val(group), path("span_diff_target_pooling.tsv", optional: true), emit: pooling
    path "*.pooled.bam", optional: true, emit: pooled_bam
    path "*.pooled.bam.bai", optional: true, emit: pooled_bai
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def caller_arg = caller_priority ? "--caller-priority '${caller_priority}'" : ''
    def pooling_dir = "${params.outdir}/03_peak_calling/08_differential/03_span/${group}"
    """
    span_diff.py \
        --jar ${omnipeaks_jar} \
        --mode ${mode} \
        --contrast '${contrast}' \
        --group '${group}' \
        --samples ${samples_manifest} \
        --peaks ${peaks_manifest} \
        ${caller_arg} \
        --chrom-sizes ${chrom_sizes} \
        --bin ${bin_size} \
        --gap ${gap} \
        --fdr ${fdr} \
        --fallback-backend ${fallback_backend} \
        --use-spikein ${use_spikein} \
        --allow-partial ${allow_partial} \
        --java-heap ${java_heap} \
        --cpus ${task.cpus} \
        --pooling-dir '${pooling_dir}' \
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o "([0-9]{1,}\\.)+[0-9]{1,}")
        java: \$(java -version 2>&1 | head -n 1 | sed -e 's/"//g')
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
        samtools: \$(samtools --version | head -n 1 | sed 's/samtools //')
        R: \$(R --version | head -n 1 | sed 's/.* //')
    END_VERSIONS
    """
}
