/*
 * Differential-only entrypoint
 */

include { DIFFERENTIAL_PEAK_CALLING } from "../subworkflows/local/differential_peak_calling"
include { CUSTOM_GETCHROMSIZES } from "../modules/nf-core/custom/getchromsizes/main"

workflow DIFFERENTIAL_ONLY {

    if (!params.differential_from_run) {
        exit 1, "--differential_from_run must be provided when using -entry DIFFERENTIAL_ONLY"
    }

    def manifest_dir = "${params.differential_from_run}/03_peak_calling/08_differential/00_manifests"
    def samples_manifest = "${manifest_dir}/differential_manifest.samples.tsv"
    def peaks_manifest = "${manifest_dir}/differential_manifest.peaks.tsv"
    def cached_gene_bed = "${manifest_dir}/differential_manifest.gene_bed.bed"

    if (!file(samples_manifest).exists()) {
        exit 1, "Missing samples manifest: ${samples_manifest}"
    }
    if (!file(peaks_manifest).exists()) {
        exit 1, "Missing peaks manifest: ${peaks_manifest}"
    }

    ch_samples_manifest = Channel.fromPath(samples_manifest, checkIfExists: true)
    ch_peaks_manifest = Channel.fromPath(peaks_manifest, checkIfExists: true)

    if (params.run_span_diff && !params.omnipeaks_jar) {
        exit 1, "--omnipeaks_jar is required when --run_span_diff is enabled."
    }

    if (!params.gene_bed && !file(cached_gene_bed).exists() && !params.gtf) {
        exit 1, "Differential annotation requires --gene_bed or --gtf (or a cached gene bed at ${cached_gene_bed})."
    }

    ch_gene_bed = Channel.empty()
    if (params.gene_bed) {
        ch_gene_bed = Channel.fromPath(params.gene_bed, checkIfExists: true)
    } else if (file(cached_gene_bed).exists()) {
        ch_gene_bed = Channel.fromPath(cached_gene_bed, checkIfExists: true)
    }

    ch_gtf = Channel.empty()
    if (params.gtf) {
        ch_gtf = Channel.from( file(params.gtf) )
    }

    ch_chrom_sizes = Channel.empty()
    def need_chrom_sizes = params.run_chipbinner || (params.run_span_diff && params.span_diff_mode != 'fallback')
    if (need_chrom_sizes) {
        if (!params.fasta) {
            exit 1, "--fasta is required to compute chrom sizes for differential-only mode when running ChIPBinner or SPAN native/auto."
        }
        ch_fasta = Channel.of([ [id: 'genome'], file(params.fasta) ])
        CUSTOM_GETCHROMSIZES ( ch_fasta )
        ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    }

    DIFFERENTIAL_PEAK_CALLING (
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        ch_chrom_sizes,
        ch_gene_bed,
        ch_gtf,
        Channel.empty(),
        ch_samples_manifest,
        ch_peaks_manifest,
        'posthoc'
    )
}
