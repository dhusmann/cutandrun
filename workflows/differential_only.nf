/*
 * Differential-only entrypoint
 */

include { DIFFERENTIAL_PEAK_CALLING } from "../subworkflows/local/differential_peak_calling"
include { CUSTOM_GETCHROMSIZES } from "../modules/nf-core/custom/getchromsizes/main"

import java.security.MessageDigest

def hashFile(path) {
    if (!path) {
        return "none"
    }
    def f = file(path)
    if (!f.exists()) {
        return "none"
    }
    def digest = MessageDigest.getInstance("SHA-256")
    f.withInputStream { input ->
        byte[] buffer = new byte[8192]
        int read
        while ((read = input.read(buffer)) > 0) {
            digest.update(buffer, 0, read)
        }
    }
    return digest.digest().encodeHex().toString().substring(0, 12)
}

def resolveCachedWindowsPath(windowsDir, binSize, genomeId, blacklistPath) {
    if (!windowsDir) {
        return null
    }
    def winPath = file(windowsDir)
    if (winPath.exists() && winPath.isFile()) {
        return winPath.toString()
    }
    if (!winPath.exists() || !winPath.isDirectory()) {
        return null
    }
    def blacklistHash = hashFile(blacklistPath)
    def meta = file("${windowsDir}/chipbinner_windows_meta.tsv")
    if (meta.exists()) {
        def lines = meta.text.readLines()
        if (lines.size() > 1) {
            def header = lines[0].split("\t")
            def values = lines[1].split("\t", -1)
            def metaRow = [:]
            header.eachWithIndex { col, idx ->
                metaRow[col] = idx < values.size() ? values[idx] : ""
            }
            def genome_ok = !genomeId || metaRow["genome_id"] == genomeId
            if (genome_ok && metaRow["bin_size"] == binSize.toString() && metaRow["blacklist_hash"] == blacklistHash) {
                def candidate = file("${windowsDir}/${metaRow['windows_path']}")
                if (candidate.exists()) {
                    return candidate.toString()
                }
            }
        }
    }
    def preferred = []
    if (genomeId) {
        preferred << "windows.${genomeId}.${binSize}.${blacklistHash}.bed"
        preferred << "windows.${genomeId}.${binSize}.bed"
    } else {
        preferred << "windows.${binSize}.${blacklistHash}.bed"
        preferred << "windows.${binSize}.bed"
    }
    for (name in preferred) {
        def candidate = file("${windowsDir}/${name}")
        if (candidate.exists()) {
            return candidate.toString()
        }
    }
    return null
}

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

    def annotation_enabled = (params.run_diffbind || params.run_chipbinner || params.run_span_diff) && !params.differential_publish_manifest_only
    if (annotation_enabled && !params.gene_bed && !file(cached_gene_bed).exists() && !params.gtf) {
        exit 1, "Differential annotation requires --gene_bed or --gtf (or a cached gene bed at ${cached_gene_bed})."
    }

    ch_gene_bed = Channel.empty()
    if (annotation_enabled && params.gene_bed) {
        ch_gene_bed = Channel.fromPath(params.gene_bed, checkIfExists: true)
    } else if (annotation_enabled && file(cached_gene_bed).exists()) {
        ch_gene_bed = Channel.fromPath(cached_gene_bed, checkIfExists: true)
    }

    ch_gtf = Channel.empty()
    if (annotation_enabled && params.gtf) {
        ch_gtf = Channel.from( file(params.gtf) )
    }

    def cached_windows = resolveCachedWindowsPath(
        params.chipbinner_windows_dir,
        params.chipbinner_bin_size,
        params.genome,
        params.blacklist
    )
    def use_cached_windows = cached_windows != null
    if (use_cached_windows) {
        params.chipbinner_windows_dir = cached_windows
    }

    ch_chrom_sizes = Channel.empty()
    def span_requires_sizes = params.run_span_diff && ['native', 'auto'].contains(params.span_diff_mode)
    def chipbinner_requires_sizes = params.run_chipbinner && !use_cached_windows
    def need_real_chrom_sizes = span_requires_sizes || chipbinner_requires_sizes
    if (need_real_chrom_sizes) {
        if (!params.fasta) {
            def reasons = []
            if (chipbinner_requires_sizes) {
                reasons << "cached ChIPBinner windows not found"
            }
            if (span_requires_sizes) {
                if (params.span_diff_mode == 'auto') {
                    reasons << "SPAN auto can select native compare and requires chrom sizes (set --span_diff_mode fallback to run without FASTA)"
                } else {
                    reasons << "SPAN native requires chrom sizes"
                }
            }
            exit 1, "--fasta is required for differential-only mode (${reasons.join('; ')})"
        }
        ch_fasta = Channel.of([ [id: 'genome'], file(params.fasta) ])
        CUSTOM_GETCHROMSIZES ( ch_fasta )
        ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    } else if (params.run_chipbinner || params.run_span_diff) {
        ch_chrom_sizes = Channel.fromPath("${projectDir}/assets/chrom_sizes_stub.sizes", checkIfExists: true)
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
