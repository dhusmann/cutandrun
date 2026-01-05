/*
 * Differential peak calling subworkflow
 */

import groovy.json.JsonOutput

include { DIFFERENTIAL_MANIFESTS } from "../../modules/local/differential_manifests"
include { DIFFERENTIAL_DESIGN } from "../../modules/local/differential_design"
include { PUBLISH_DIFFERENTIAL_MANIFESTS } from "../../modules/local/publish_differential_manifests"
include { RECORDS_TO_TSV as RECORDS_TO_TSV_DIFFBIND } from "../../modules/local/records_to_tsv"
include { RECORDS_TO_TSV as RECORDS_TO_TSV_CHIPBINNER } from "../../modules/local/records_to_tsv"
include { RECORDS_TO_TSV as RECORDS_TO_TSV_SPAN } from "../../modules/local/records_to_tsv"
include { DIFFBIND_RUN } from "../../modules/local/diffbind_run"
include { CHIPBINNER_RUN } from "../../modules/local/chipbinner_run"
include { SPAN_DIFF_RUN } from "../../modules/local/span_diff_run"
include { ANNOTATE_REGIONS } from "../../modules/local/annotate_regions"
include { DIFFERENTIAL_SUMMARY } from "../../modules/local/differential_summary"
include { MULTIQC_DIFFERENTIAL } from "../../modules/local/multiqc_differential"

workflow DIFFERENTIAL_PEAK_CALLING {
    take:
    ch_bam
    ch_bai
    ch_bigwig
    ch_peaks
    ch_scale_factors
    ch_chrom_sizes
    ch_gene_bed
    ch_blacklist
    ch_samples_manifest_in
    ch_peaks_manifest_in
    mode

    main:
    ch_versions = Channel.empty()
    ch_diff_summary_mqc = Channel.empty()
    ch_diff_design_mqc = Channel.empty()
    ch_annotation_requests = Channel.empty()
    ch_diffbind_summary = Channel.empty()
    ch_chipbinner_summary = Channel.empty()
    ch_span_summary = Channel.empty()
    ch_summary_files = Channel.empty()
    ch_samples_manifest = Channel.empty()
    ch_peaks_manifest = Channel.empty()
    summary_stub = file("$projectDir/assets/differential_summary_stub.tsv")

    def diff_enabled = params.run_diffbind || params.run_chipbinner || params.run_span_diff || params.differential_publish_manifest_only
    if (diff_enabled) {
        def use_manifests = (mode == 'posthoc')
        def diff_outdir = "${params.outdir}/03_peak_calling/08_differential"

        def bam_subdir = params.run_remove_linear_dups ? 'linear_dedup' : (params.run_remove_dups ? 'dedup' : (params.run_mark_dups ? 'markdup' : ''))
        def bam_dir = bam_subdir ? "${params.outdir}/02_alignment/${params.aligner}/target/${bam_subdir}" : "${params.outdir}/02_alignment/${params.aligner}/target"

        def peaks_dir = { caller -> "${params.outdir}/03_peak_calling/04_called_peaks/${caller}" }
        def bigwig_dir = "${params.outdir}/03_peak_calling/03_bed_to_bigwig"

        def peaks_format = { file ->
            def name = file.getName().toLowerCase()
            if (name.endsWith('.narrowpeak')) return 'narrowPeak'
            if (name.endsWith('.broadpeak')) return 'broadPeak'
            if (name.endsWith('.peak')) return 'peak'
            if (name.endsWith('.bed')) return 'bed'
            return 'bed'
        }

        if (!use_manifests) {
            ch_bam_target = ch_bam.filter { it[0].is_control == false }
            ch_bai_target = ch_bai.filter { it[0].is_control == false }
            ch_bigwig_target = ch_bigwig.filter { it[0].is_control == false }
            ch_scale_target = ch_scale_factors.filter { it[0].is_control == false }

            ch_scale_map = ch_scale_target
                .map { meta, scale -> [meta.id, scale] }
                .collect()
                .map { items ->
                    def out = [:]
                    items.each { entry ->
                        if (entry instanceof List && entry.size() > 1) {
                            out[entry[0]] = entry[1]
                        }
                    }
                    out
                }

            ch_bigwig_map = ch_bigwig_target
                .map { meta, bw -> [meta.id, bw] }
                .collect()
                .map { items ->
                    def out = [:]
                    items.each { entry ->
                        if (entry instanceof List && entry.size() > 1) {
                            out[entry[0]] = entry[1]
                        }
                    }
                    out
                }

            ch_samples_raw = ch_bam_target
                .map { meta, bam -> [meta.id, meta, bam] }
                .join(ch_bai_target.map { meta, bai -> [meta.id, bai] })
                .map { id, meta, bam, bai -> [meta, bam, bai] }
                .combine(ch_bigwig_map)
                .combine(ch_scale_map)
                .map { meta, bam, bai, bw_map, scale_map ->
                    def bigwig = bw_map.get(meta.id) ?: 'NA'
                    def scale = scale_map.get(meta.id) ?: 'NA'
                    def bam_path = "${bam_dir}/${bam.getName()}"
                    def bai_path = "${bam_dir}/${bai.getName()}"
                    def bigwig_path = bigwig == 'NA' ? 'NA' : "${bigwig_dir}/${bigwig.getName()}"
                    [
                        sample_id: meta.id,
                        group: meta.group,
                        condition: meta.condition,
                        replicate: meta.replicate,
                        final_bam: bam_path,
                        final_bai: bai_path,
                        spikein_scale_factor: scale,
                        bigwig_path: bigwig_path
                    ]
                }
                .map { row ->
                    "${row.sample_id}\t${row.group}\t${row.condition}\t${row.replicate}\t${row.final_bam}\t${row.final_bai}\t${row.spikein_scale_factor}\t${row.bigwig_path}"
                }

            def samples_header = "sample_id\tgroup\tcondition\treplicate\tfinal_bam\tfinal_bai\tspikein_scale_factor\tbigwig_path"
            ch_samples_raw = ch_samples_raw
                .collect()
                .map { rows -> ([samples_header] + rows).join('\n') + '\n' }
                .collectFile(name: 'differential_samples_raw.tsv', newLine: false)

            ch_peaks_raw = ch_peaks
                .map { meta, peak_file ->
                    def peak_path = "${peaks_dir(meta.caller)}/${peak_file.getName()}"
                    def fmt = peaks_format(peak_file)
                    "${meta.id}\t${meta.group}\t${meta.condition}\t${meta.replicate}\t${meta.caller}\t${peak_path}\t${fmt}"
                }

            def peaks_header = "sample_id\tgroup\tcondition\treplicate\tcaller\tpeaks_path\tpeaks_format"
            ch_peaks_raw = ch_peaks_raw
                .collect()
                .map { rows -> ([peaks_header] + rows).join('\n') + '\n' }
                .collectFile(name: 'differential_peaks_raw.tsv', newLine: false)

            DIFFERENTIAL_MANIFESTS (
                ch_samples_raw,
                ch_peaks_raw,
                params.chipbinner_ms_coeffs ? file(params.chipbinner_ms_coeffs) : file("$projectDir/assets/ms_coeffs_stub.tsv"),
                params.normalisation_mode
            )
            ch_versions = ch_versions.mix(DIFFERENTIAL_MANIFESTS.out.versions)
            ch_samples_manifest = DIFFERENTIAL_MANIFESTS.out.samples
            ch_peaks_manifest = DIFFERENTIAL_MANIFESTS.out.peaks
        } else {
            ch_samples_manifest = ch_samples_manifest_in
            ch_peaks_manifest = ch_peaks_manifest_in
        }

        DIFFERENTIAL_DESIGN (
            ch_samples_manifest,
            ch_peaks_manifest,
            params.differential_contrast,
            params.differential_min_replicates,
            params.differential_allow_partial,
            params.run_diffbind,
            params.run_chipbinner,
            params.run_span_diff,
            params.differential_groups ?: '',
            params.differential_callers ?: ''
        )
        ch_versions = ch_versions.mix(DIFFERENTIAL_DESIGN.out.versions)

        PUBLISH_DIFFERENTIAL_MANIFESTS (
            ch_samples_manifest,
            ch_peaks_manifest
        )
        ch_versions = ch_versions.mix(PUBLISH_DIFFERENTIAL_MANIFESTS.out.versions)

        ch_samples_rows = ch_samples_manifest.splitCsv(header: true, sep: '\t')
        ch_peaks_rows = ch_peaks_manifest.splitCsv(header: true, sep: '\t')
        ch_design_rows = DIFFERENTIAL_DESIGN.out.design.splitCsv(header: true, sep: '\t')

        ch_gene_bed_single = ch_gene_bed.first()

        def contrast_labels = params.differential_contrast.split(',').collect { it.trim() }
        ch_samples_rows_contrast = ch_samples_rows.filter { row -> contrast_labels.contains(row.condition) }

        def diff_use_spikein = params.differential_use_spikein ? params.differential_use_spikein.toString().toLowerCase() : 'auto'
        def use_spikein = (diff_use_spikein == 'true') || (diff_use_spikein == 'auto' && params.normalisation_mode == 'Spikein')

        if (params.run_diffbind && !params.differential_publish_manifest_only) {
            ch_diffbind_design = ch_design_rows
                .filter { row -> row.caller != 'NA' && row.status == 'RUN' && row.eligible_diffbind == 'true' }
                .map { row -> [ [row.group, row.caller], row ] }

            ch_sample_peak_records = ch_samples_rows_contrast
                .map { row -> [row.sample_id, row] }
                .join(ch_peaks_rows.map { row -> [row.sample_id, row] })
                .map { sample_id, srow, prow ->
                    def record = [
                        group: srow.group,
                        caller: prow.caller,
                        sample_id: srow.sample_id,
                        condition: srow.condition,
                        replicate: srow.replicate,
                        bam: srow.final_bam,
                        peaks: prow.peaks_path,
                        spikein_scale_factor: srow.spikein_scale_factor
                    ]
                    [ [record.group, record.caller], record ]
                }
                .groupTuple(by: [0])
                .map { key, recs -> [ key[0], key[1], recs ] }

            ch_diffbind_records = ch_sample_peak_records
                .map { group, caller, recs -> [ [group, caller], recs ] }
                .join(ch_diffbind_design)
                .map { key, recs, row -> [ key[0], key[1], recs ] }

            ch_diffbind_records
                .map { group, caller, recs -> [ group, caller, JsonOutput.toJson(recs) ] }
                .set { ch_diffbind_records_json }

            RECORDS_TO_TSV_DIFFBIND (
                ch_diffbind_records_json.map { it[0] },
                ch_diffbind_records_json.map { it[1] },
                ch_diffbind_records_json.map { it[2] },
                'sample_id\tgroup\tcondition\treplicate\tbam\tpeaks\tcaller\tspikein_scale_factor',
                ch_diffbind_records_json.map { group, caller, json -> "diffbind_${caller}_${group}.tsv" }
            )
            ch_versions = ch_versions.mix(RECORDS_TO_TSV_DIFFBIND.out.versions)

            DIFFBIND_RUN (
                RECORDS_TO_TSV_DIFFBIND.out.tsv,
                params.differential_contrast,
                use_spikein,
                params.diffbind_fdr,
                params.diffbind_lfc,
                params.diffbind_min_overlap,
                params.diffbind_backend,
                params.diffbind_recenter_peaks,
                params.diffbind_summits,
                params.diffbind_norm_method,
                params.diffbind_extra_params ? file(params.diffbind_extra_params).toString() : '',
                params.export_diffbind_sheets
            )
            ch_versions = ch_versions.mix(DIFFBIND_RUN.out.versions)

            ch_diffbind_annot = DIFFBIND_RUN.out.results
                .combine(ch_gene_bed_single)
                .map { tuple, gene_bed ->
                    def group = tuple[0]
                    def caller = tuple[1]
                    def regions = tuple[2]
                    [ 'diffbind', group, caller, regions, gene_bed, "01_diffbind/${caller}/${group}/diffbind.results.annotated.tsv" ]
                }

            ch_annotation_requests = ch_annotation_requests.mix(ch_diffbind_annot)
            ch_diffbind_summary = DIFFBIND_RUN.out.summary.map { it[2] }
        }

        if (params.run_chipbinner && !params.differential_publish_manifest_only) {
            ch_chip_design = ch_design_rows
                .filter { row -> row.caller == 'NA' && row.status == 'RUN' && row.eligible_chipbinner == 'true' }
                .map { row -> [row.group, row] }

            ch_chip_records = ch_samples_rows_contrast
                .map { row -> [row.group, row] }
                .groupTuple(by: [0])
                .join(ch_chip_design)
                .map { group, rows, design -> [ group, rows ] }

            ch_chip_records
                .map { group, rows -> [ group, JsonOutput.toJson(rows) ] }
                .set { ch_chip_records_json }

            RECORDS_TO_TSV_CHIPBINNER (
                ch_chip_records_json.map { it[0] },
                ch_chip_records_json.map { 'NA' },
                ch_chip_records_json.map { it[1] },
                'sample_id\tgroup\tcondition\treplicate\tfinal_bam\tfinal_bai\tbigwig_path\tspikein_scale_factor\tms_coeff',
                ch_chip_records_json.map { group, json -> "chipbinner_${group}.tsv" }
            )
            ch_versions = ch_versions.mix(RECORDS_TO_TSV_CHIPBINNER.out.versions)

            ch_chip_records_file = RECORDS_TO_TSV_CHIPBINNER.out.tsv
                .map { group, caller, records_file -> [group, records_file] }

            ch_chrom_sizes_single = ch_chrom_sizes.collect().map { it[0] }

            ch_chip_inputs = ch_chip_records_file
                .combine(ch_chrom_sizes_single)
                .map { record, chrom_sizes -> [ record[0], record[1], chrom_sizes ] }

            CHIPBINNER_RUN (
                ch_chip_inputs,
                params.differential_contrast,
                params.chipbinner_bin_size,
                params.chipbinner_windows_dir ? file(params.chipbinner_windows_dir).toString() : '',
                params.blacklist ? file(params.blacklist).toString() : '',
                params.chipbinner_use_input,
                use_spikein,
                params.chipbinner_pseudocount,
                params.chipbinner_hdbscan_grid_minpts,
                params.chipbinner_hdbscan_grid_minsamps,
                params.chipbinner_fdr,
                params.chipbinner_lfc,
                params.chipbinner_bootstrap,
                params.chipbinner_k_value,
                params.chipbinner_functional_db ? file(params.chipbinner_functional_db).toString() : '',
                params.differential_allow_partial
            )
            ch_versions = ch_versions.mix(CHIPBINNER_RUN.out.versions)

            ch_chip_annot = CHIPBINNER_RUN.out.differential
                .combine(ch_gene_bed_single)
                .map { tuple, gene_bed ->
                    def group = tuple[0]
                    def regions = tuple[1]
                    [ 'chipbinner', group, 'NA', regions, gene_bed, "02_chipbinner/${group}/chipbinner.differential.annotated.tsv" ]
                }

            ch_annotation_requests = ch_annotation_requests.mix(ch_chip_annot)
            ch_chipbinner_summary = CHIPBINNER_RUN.out.summary.map { it[1] }
        }

        if (params.run_span_diff && !params.differential_publish_manifest_only) {
            ch_span_design = ch_design_rows
                .filter { row -> row.caller == 'NA' && row.status == 'RUN' && row.eligible_span == 'true' }
                .map { row -> [row.group, row] }

            ch_span_records = ch_samples_rows_contrast
                .map { row -> [row.group, row] }
                .groupTuple(by: [0])
                .join(ch_span_design)
                .map { group, rows, design -> [ group, rows ] }

            ch_span_records
                .map { group, rows -> [ group, JsonOutput.toJson(rows) ] }
                .set { ch_span_records_json }

            RECORDS_TO_TSV_SPAN (
                ch_span_records_json.map { it[0] },
                ch_span_records_json.map { 'NA' },
                ch_span_records_json.map { it[1] },
                'sample_id\tgroup\tcondition\treplicate\tfinal_bam',
                ch_span_records_json.map { group, json -> "span_${group}.tsv" }
            )
            ch_versions = ch_versions.mix(RECORDS_TO_TSV_SPAN.out.versions)

            ch_span_records_file = RECORDS_TO_TSV_SPAN.out.tsv
                .map { group, caller, records_file -> [group, records_file] }

            def span_caller_priority = (params.callers ?: [])
                .collect { it.toString().toLowerCase() }
                .findAll { it.startsWith('span') || it.startsWith('omnipeak') }

            ch_span_peaks = ch_peaks_rows
                .filter { row ->
                    def caller = row.caller?.toString()?.toLowerCase()
                    caller && (caller.startsWith('span') || caller.startsWith('omnipeak'))
                }
                .map { row -> [row.group, row.caller?.toString()?.toLowerCase(), row.peaks_path] }
                .groupTuple(by: [0])
                .map { group, entries ->
                    def peaks_by_caller = entries.groupBy { it[1] }
                    def chosen = span_caller_priority.find { peaks_by_caller.containsKey(it) } ?: peaks_by_caller.keySet().sort()[0]
                    if (peaks_by_caller.size() > 1) {
                        log.warn "Multiple SPAN callers for group ${group}; using '${chosen}' for differential peaks."
                    }
                    [group, peaks_by_caller[chosen][0][2]]
                }

            ch_span_inputs = ch_span_records_file
                .join(ch_span_peaks)
                .map { group, records_file, peaks_path -> [ group, records_file, peaks_path ] }

            SPAN_DIFF_RUN (
                ch_span_inputs,
                params.differential_contrast,
                params.span_diff_mode,
                params.span_diff_fdr,
                params.span_diff_gap,
                params.span_diff_bin,
                params.span_fallback_backend,
                params.omnipeaks_jar ? file(params.omnipeaks_jar) : null,
                params.span_diff_java_heap
            )
            ch_versions = ch_versions.mix(SPAN_DIFF_RUN.out.versions)

            ch_span_annot = SPAN_DIFF_RUN.out.differential
                .combine(ch_gene_bed_single)
                .map { tuple, gene_bed ->
                    def group = tuple[0]
                    def regions = tuple[1]
                    [ 'span', group, 'NA', regions, gene_bed, "03_span/${group}/span.differential.annotated.tsv" ]
                }

            ch_annotation_requests = ch_annotation_requests.mix(ch_span_annot)
            ch_span_summary = SPAN_DIFF_RUN.out.summary.map { it[1] }
        }

        ANNOTATE_REGIONS (
            ch_annotation_requests
        )
        ch_versions = ch_versions.mix(ANNOTATE_REGIONS.out.versions)

        ch_summary_files = ch_diffbind_summary.mix(ch_chipbinner_summary).mix(ch_span_summary)
        DIFFERENTIAL_SUMMARY (
            DIFFERENTIAL_DESIGN.out.design,
            ch_summary_files.collect().map { it + summary_stub }.ifEmpty([summary_stub]),
            params.differential_publish_manifest_only,
            file("$projectDir/assets/multiqc/differential_summary_header.txt"),
            file("$projectDir/assets/multiqc/differential_design_header.txt")
        )
        ch_versions = ch_versions.mix(DIFFERENTIAL_SUMMARY.out.versions)
        ch_diff_summary_mqc = DIFFERENTIAL_SUMMARY.out.summary
        ch_diff_design_mqc = DIFFERENTIAL_SUMMARY.out.design

        if (params.differential_multiqc_report) {
            MULTIQC_DIFFERENTIAL (
                file("$projectDir/assets/multiqc_differential_config.yml"),
                diff_outdir
            )
            ch_versions = ch_versions.mix(MULTIQC_DIFFERENTIAL.out.versions)
        }
    }

    emit:
    versions = ch_versions
    summary_mqc = ch_diff_summary_mqc
    design_mqc = ch_diff_design_mqc
}
