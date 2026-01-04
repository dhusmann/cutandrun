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

    ch_samples_manifest = Channel.fromPath(samples_manifest, checkIfExists: true)
    ch_peaks_manifest = Channel.fromPath(peaks_manifest, checkIfExists: true)

    if (!params.gene_bed) {
        exit 1, "--gene_bed (or --genome with a bed12 attribute) is required for differential analysis."
    }
    if (params.run_span_diff && !params.omnipeaks_jar) {
        exit 1, "--omnipeaks_jar is required when --run_span_diff is enabled."
    }

    ch_gene_bed = Channel.fromPath(params.gene_bed, checkIfExists: true)

    ch_chrom_sizes = Channel.empty()
    if (params.run_chipbinner) {
        if (!params.fasta) {
            exit 1, "--fasta is required to compute chrom sizes for ChIPBinner in differential-only mode."
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
        Channel.empty(),
        ch_samples_manifest,
        ch_peaks_manifest,
        'posthoc'
    )
}
