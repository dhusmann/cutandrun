/*
 * Differential-only entrypoint (posthoc)
 */

import groovy.json.JsonSlurper

include { DIFFERENTIAL_PEAK_CALLING } from '../subworkflows/local/differential_peak_calling'

workflow DIFFERENTIAL_ONLY {
    main:
    if (!params.differential_from_run) {
        exit 1, "--differential_from_run must be provided for DIFFERENTIAL_ONLY"
    }

    def run_dir = params.differential_from_run
    def samples_path = file("${run_dir}/03_peak_calling/06_differential/00_manifests/differential_manifest.samples.tsv")
    def peaks_path = file("${run_dir}/03_peak_calling/06_differential/00_manifests/differential_manifest.peaks.tsv")
    def run_meta_path = file("${run_dir}/03_peak_calling/06_differential/00_manifests/differential_manifest.run_meta.json")
    def chrom_sizes_path = file("${run_dir}/03_peak_calling/06_differential/00_manifests/chrom_sizes.sizes")

    if (!samples_path.exists()) {
        exit 1, "Missing differential manifest samples.tsv under ${run_dir}"
    }
    if (!peaks_path.exists()) {
        exit 1, "Missing differential manifest peaks.tsv under ${run_dir}"
    }
    if (run_meta_path.exists() && !params.differential_contrast) {
        def meta = new JsonSlurper().parse(run_meta_path)
        if (meta?.default_contrast) {
            params.differential_contrast = meta.default_contrast
        }
    }

    def callers_list = []
    if (run_meta_path.exists()) {
        def meta = new JsonSlurper().parse(run_meta_path)
        callers_list = meta?.callers ?: []
    }

    ch_samples = Channel.fromPath(samples_path).splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [
                id: row.sample_id,
                sample_id: row.sample_id,
                group: row.group,
                condition: row.condition,
                replicate: row.replicate.toInteger(),
                is_control: false
            ]
            def bam = file(row.final_bam)
            def bai = file(row.final_bai)
            [meta, bam, bai, row.spikein_scale_factor, row.bigwig]
        }

    ch_bam_bai = ch_samples.map { meta, bam, bai, scale, bigwig -> [meta, bam, bai] }
    ch_spikein_scale = ch_samples
        .filter { meta, bam, bai, scale, bigwig -> scale && scale != 'NA' }
        .map { meta, bam, bai, scale, bigwig -> [meta, scale] }

    ch_bigwig = ch_samples
        .filter { meta, bam, bai, scale, bigwig -> bigwig && bigwig != 'NA' }
        .map { meta, bam, bai, scale, bigwig -> [meta, file(bigwig)] }

    ch_samples_map = ch_samples
        .map { meta, bam, bai, scale, bigwig -> [meta.sample_id, meta] }
        .toList()
        .map { list -> list.collectEntries { [(it[0]): it[1]] } }

    ch_peaks = Channel.fromPath(peaks_path).splitCsv(header: true, sep: '\t')
        .combine(ch_samples_map)
        .map { row, samples_map ->
            def meta = samples_map[row.sample_id]
            if (!meta) {
                return null
            }
            def meta_out = meta + [caller: row.caller]
            [meta_out, file(row.peaks_path)]
        }
        .filter { it != null }

    if (!callers_list) {
        callers_list = peaks_path.text.readLines()
            .drop(1)
            .collect { it.split('\\t')[1] }
            .unique()
    }

    ch_chrom_sizes = chrom_sizes_path.exists() ? Channel.value(chrom_sizes_path) : Channel.empty()

    DIFFERENTIAL_PEAK_CALLING(
        ch_bam_bai,
        ch_peaks,
        ch_bigwig,
        ch_spikein_scale,
        ch_chrom_sizes,
        callers_list
    )
    ch_diff_summary = DIFFERENTIAL_PEAK_CALLING.out.summary
    ch_diff_skipped = DIFFERENTIAL_PEAK_CALLING.out.skipped
    ch_diff_versions = DIFFERENTIAL_PEAK_CALLING.out.versions
    if (!params.run_multiqc) {
        ch_diff_summary.subscribe { }
        ch_diff_skipped.subscribe { }
        ch_diff_versions.subscribe { }
    }
}
