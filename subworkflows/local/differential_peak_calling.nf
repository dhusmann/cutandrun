/*
 * Differential peak calling / enrichment
 */

import groovy.json.JsonOutput
import groovy.json.JsonSlurper

include { WRITE_DIFFERENTIAL_MANIFESTS } from '../../modules/local/write_differential_manifests'
include { MAKE_DIFFBIND_SAMPLESHEET } from '../../modules/local/diffbind_samplesheet'
include { MAKE_DIFFBIND_SAMPLESHEET as MAKE_SPAN_SAMPLESHEET } from '../../modules/local/diffbind_samplesheet'
include { DIFFBIND_RUN } from '../../modules/local/diffbind_run'
include { CHIPBINNER_BINS } from '../../modules/local/chipbinner_bins'
include { CHIPBINNER_COUNTS } from '../../modules/local/chipbinner_counts'
include { CHIPBINNER_HDBSCAN_GRID } from '../../modules/local/chipbinner_hdbscan_grid'
include { CHIPBINNER_ROTS } from '../../modules/local/chipbinner_rots'
include { SPAN_CAPABILITY_PROBE } from '../../modules/local/span_capability_probe'
include { SPAN_COMPARE } from '../../modules/local/span_compare'
include { SPAN_FALLBACK_DIFF } from '../../modules/local/span_fallback_diff'
include { ANNOTATE_REGIONS } from '../../modules/local/annotate_regions'
include { DIFFERENTIAL_OVERLAP } from '../../modules/local/differential_overlap'
include { DIFFERENTIAL_SUMMARY_MERGE } from '../../modules/local/differential_summary'

workflow DIFFERENTIAL_PEAK_CALLING {
    take:
    ch_bam_bai           // tuple [meta, bam, bai]
    ch_peaks             // tuple [meta, peaks]
    ch_bigwig            // tuple [meta, bigwig]
    ch_spikein_scale     // tuple [meta, scale]
    ch_chrom_sizes       // path
    callers              // list of callers

    main:
    ch_versions = Channel.empty()
    ch_skipped_records = Channel.empty()
    ch_summary_files = Channel.empty()
    ch_diffbind_significant = Channel.empty()
    ch_chipbinner_significant = Channel.empty()
    ch_span_significant = Channel.empty()

    ch_chrom_sizes_single = ch_chrom_sizes
        .map { it instanceof List ? it[0] : it }
        .ifEmpty('')
    ch_gene_bed = Channel.value(params.gene_bed ? file(params.gene_bed) : '')
    ch_gtf = Channel.value(params.gtf ? file(params.gtf) : '')

    def contrast_parts = params.differential_contrast?.split(',')?.collect { it.trim() }?.findAll { it }
    if (!contrast_parts || contrast_parts.size() != 2) {
        exit 1, "Invalid --differential_contrast value. Expected 'TREATED,CONTROL'."
    }
    def treated = contrast_parts[0]
    def control = contrast_parts[1]

    def use_spikein_param = params.differential_use_spikein?.toString()?.toLowerCase()
    ch_use_spikein = Channel.value('false')
    if (use_spikein_param == 'true') {
        ch_use_spikein = Channel.value('true')
    } else if (use_spikein_param == 'auto') {
        ch_use_spikein = ch_spikein_scale
            .map { meta, scale -> scale }
            .filter { it != null && it.toString() != 'NA' }
            .count()
            .map { it > 0 ? 'true' : 'false' }
    }

    ch_bam_bai_target = ch_bam_bai.filter { meta, bam, bai -> meta.is_control == false }
    ch_bigwig_target = ch_bigwig.filter { meta, bw -> meta.is_control == false }
    ch_scale_target = ch_spikein_scale.filter { meta, scale -> meta.is_control == false }

    ch_bigwig_map = ch_bigwig_target
        .map { meta, bw -> [meta.sample_id ?: meta.id, bw] }
        .toList()
        .map { list -> list.collectEntries { [(it[0]): it[1]] } }

    ch_scale_map = ch_scale_target
        .map { meta, scale -> [meta.sample_id ?: meta.id, scale] }
        .toList()
        .map { list -> list.collectEntries { [(it[0]): it[1]] } }

    ch_samples = ch_bam_bai_target
        .combine(ch_bigwig_map)
        .combine(ch_scale_map)
        .map { meta, bam, bai, bw_map, scale_map ->
            def sample_id = meta.sample_id ?: meta.id
            def scale = scale_map.get(sample_id, 'NA')
            def bigwig = bw_map.get(sample_id, 'NA')
            [
                sample_id: sample_id,
                group: meta.group,
                condition: meta.condition,
                replicate: meta.replicate,
                bam: bam,
                bai: bai,
                bam_path: bam.toString(),
                bai_path: bai.toString(),
                normalisation_mode: params.normalisation_mode,
                spikein_scale_factor: scale,
                bigwig_path: bigwig == null ? 'NA' : bigwig.toString(),
                meta: meta
            ]
        }

    // Write manifests in integrated mode
    if (!params.differential_from_run) {
        ch_samples_manifest = ch_samples
            .map { rec ->
                [
                    sample_id: rec.sample_id,
                    group: rec.group,
                    condition: rec.condition,
                    replicate: rec.replicate,
                    final_bam: rec.bam_path,
                    final_bai: rec.bai_path,
                    normalisation_mode: rec.normalisation_mode,
                    spikein_scale_factor: rec.spikein_scale_factor ?: 'NA',
                    bigwig: rec.bigwig_path ?: 'NA'
                ]
            }
            .toList()
            .map { list -> JsonOutput.toJson(list) }

        ch_peaks_manifest = ch_peaks
            .filter { meta, peaks -> meta.is_control == false }
            .map { meta, peaks ->
                def sample_id = meta.sample_id ?: meta.id
                def caller_role = callers && callers[0] == meta.caller ? 'primary' : 'secondary'
                def peaks_format = 'other'
                def name = peaks.getName()
                if (name.endsWith('.narrowPeak')) {
                    peaks_format = 'narrowPeak'
                } else if (name.endsWith('.broadPeak')) {
                    peaks_format = 'broadPeak'
                } else if (name.endsWith('.bed')) {
                    peaks_format = 'bed'
                }
                [
                    sample_id: sample_id,
                    caller: meta.caller,
                    peaks_path: peaks.toString(),
                    peaks_format: peaks_format,
                    caller_role: caller_role,
                    group: meta.group,
                    condition: meta.condition,
                    replicate: meta.replicate
                ]
            }
            .toList()
            .map { list -> JsonOutput.toJson(list) }

        def run_meta = [
            pipeline: workflow.manifest.name,
            pipeline_version: workflow.manifest.version,
            pipeline_revision: workflow.revision ?: 'NA',
            genome: params.genome,
            callers: callers,
            default_contrast: params.differential_contrast,
            chrom_sizes: "${params.outdir}/03_peak_calling/06_differential/00_manifests/chrom_sizes.sizes"
        ]

        ch_run_meta = Channel.value(JsonOutput.toJson(run_meta))
        ch_outdir = Channel.value(params.outdir)

        WRITE_DIFFERENTIAL_MANIFESTS(
            ch_samples_manifest,
            ch_peaks_manifest,
            ch_run_meta,
            ch_outdir,
            ch_chrom_sizes_single
        )
        ch_versions = ch_versions.mix(WRITE_DIFFERENTIAL_MANIFESTS.out.versions)
    }

    // Group-level validation
    ch_group_samples = ch_samples
        .map { rec -> [rec.group, rec] }
        .groupTuple()

    ch_group_eval = ch_group_samples
        .map { group, records ->
            def treated_records = records.findAll { it.condition == treated }
            def control_records = records.findAll { it.condition == control }
            def conditions = records.collect { it.condition }.unique()
            def other_conditions = conditions.findAll { !(it in [treated, control]) }
            def reasons = []
            if (!treated_records) {
                reasons << "missing_treated"
            }
            if (!control_records) {
                reasons << "missing_control"
            }
            if (treated_records.size() < params.differential_min_replicates) {
                reasons << "treated_replicates_lt_${params.differential_min_replicates}"
            }
            if (control_records.size() < params.differential_min_replicates) {
                reasons << "control_replicates_lt_${params.differential_min_replicates}"
            }
            if (other_conditions) {
                if (params.differential_strict) {
                    exit 1, "Group '${group}' has additional conditions (${other_conditions.join(',')}) but strict mode is enabled."
                }
                log.warn "Group '${group}' has additional conditions (${other_conditions.join(',')}); ignoring for differential analysis."
            }
            def valid = reasons.isEmpty()
            if (!valid && params.differential_strict) {
                exit 1, "Group '${group}' failed differential validation: ${reasons.join(';')}"
            }
            def filtered = records.findAll { it.condition in [treated, control] }
            [group: group, records: filtered, valid: valid, reason: reasons.join(';'), details: other_conditions ? "ignored_conditions=${other_conditions.join(',')}" : ""]
        }

    ch_group_eval.branch { valid: it.valid; invalid: !it.valid }.set { ch_group_split }
    ch_valid_groups = ch_group_split.valid
    ch_invalid_groups = ch_group_split.invalid

    ch_skipped_records = ch_skipped_records.mix(
        ch_invalid_groups.map { rec -> [method: 'all', group: rec.group, caller: 'NA', reason: rec.reason, details: rec.details] }
    )

    // DiffBind per caller x group
    if (params.run_diffbind) {
        ch_peaks_map = ch_peaks
            .filter { meta, peaks -> meta.is_control == false }
            .map { meta, peaks -> ["${meta.group}__${meta.caller}", [meta, peaks]] }
            .groupTuple()
            .toList()
            .map { list -> list.collectEntries { [(it[0]): it[1]] } }

        ch_diffbind_inputs = ch_valid_groups
            .flatMap { rec ->
                callers.collect { caller -> [group: rec.group, caller: caller, records: rec.records] }
            }
            .combine(ch_peaks_map)
            .map { rec, peaks_map ->
                def key = "${rec.group}__${rec.caller}"
                def peaks_list = peaks_map.get(key, [])
                def peaks_by_sample = peaks_list.collectEntries { [(it[0].sample_id ?: it[0].id): it[1]] }
                def missing = rec.records.findAll { !peaks_by_sample.containsKey(it.sample_id) }.collect { it.sample_id }
                def valid = missing.isEmpty()
                def samples = rec.records.collect { row ->
                    [
                        sample_id: row.sample_id,
                        group: row.group,
                        condition: row.condition,
                        replicate: row.replicate,
                        bam: row.bam_path,
                        peaks: peaks_by_sample.get(row.sample_id)?.toString(),
                        caller: rec.caller,
                        spikein_scale_factor: row.spikein_scale_factor ?: 'NA'
                    ]
                }
                [group: rec.group, caller: rec.caller, samples: samples, valid: valid, missing: missing]
            }

        ch_diffbind_inputs.branch { valid: it.valid; invalid: !it.valid }.set { ch_diffbind_split }
        ch_diffbind_valid = ch_diffbind_split.valid
        ch_diffbind_invalid = ch_diffbind_split.invalid

        ch_skipped_records = ch_skipped_records.mix(
            ch_diffbind_invalid.map { rec -> [method: 'diffbind', group: rec.group, caller: rec.caller, reason: 'missing_peaks', details: rec.missing.join(',')] }
        )

        ch_diffbind_valid
            .map { rec -> [JsonOutput.toJson(rec.samples), rec.group, rec.caller] }
            .set { ch_diffbind_samples }

        MAKE_DIFFBIND_SAMPLESHEET(
            ch_diffbind_samples.map { samples_json, group, caller -> samples_json },
            ch_diffbind_samples.map { samples_json, group, caller -> group },
            ch_diffbind_samples.map { samples_json, group, caller -> caller }
        )
        ch_versions = ch_versions.mix(MAKE_DIFFBIND_SAMPLESHEET.out.versions)

        DIFFBIND_RUN(
            MAKE_DIFFBIND_SAMPLESHEET.out.samplesheet.map { group, caller, sheet -> sheet },
            Channel.value(treated),
            Channel.value(control),
            MAKE_DIFFBIND_SAMPLESHEET.out.samplesheet.map { group, caller, sheet -> group },
            MAKE_DIFFBIND_SAMPLESHEET.out.samplesheet.map { group, caller, sheet -> caller },
            ch_use_spikein,
            Channel.value(params.diffbind_fdr),
            Channel.value(params.diffbind_lfc),
            Channel.value(params.diffbind_min_overlap),
            Channel.value(params.diffbind_summits),
            Channel.value(params.diffbind_backend),
            Channel.value(params.diffbind_extra_params ? file(params.diffbind_extra_params) : ''),
            Channel.value('diffbind')
        )
        ch_versions = ch_versions.mix(DIFFBIND_RUN.out.versions)
        ch_summary_files = ch_summary_files.mix(DIFFBIND_RUN.out.summary.map { group, caller, file -> file })
        ch_diffbind_significant = DIFFBIND_RUN.out.significant

        if (params.differential_annotate) {
            ch_diffbind_annotate = DIFFBIND_RUN.out.results.map { group, caller, file ->
                def meta = [publish_dir: "${params.outdir}/03_peak_calling/06_differential/01_diffbind/${caller}/${group}"]
                [meta, file]
            }
            ANNOTATE_REGIONS(
                ch_diffbind_annotate,
                ch_gene_bed,
                ch_gtf
            )
            ch_versions = ch_versions.mix(ANNOTATE_REGIONS.out.versions)
        }
    }

    // ChIPBinner per group
    if (params.run_chipbinner) {
        CHIPBINNER_BINS(
            ch_chrom_sizes_single,
            Channel.value(params.chipbinner_bin_size),
            ch_valid_groups.map { rec -> rec.group },
            Channel.value(params.chipbinner_windows_dir ? file(params.chipbinner_windows_dir) : ''),
            Channel.value(params.chipbinner_blacklist ? file(params.chipbinner_blacklist) : '')
        )
        ch_versions = ch_versions.mix(CHIPBINNER_BINS.out.versions)

        ch_chipbinner_samples = ch_valid_groups.map { rec ->
            def ordered = rec.records.sort { a, b ->
                a.condition <=> b.condition ?: a.replicate <=> b.replicate ?: a.sample_id <=> b.sample_id
            }
            def samples = ordered.collect { row ->
                [
                    sample_id: row.sample_id,
                    condition: row.condition,
                    bam: row.bam_path,
                    spikein_scale_factor: row.spikein_scale_factor ?: 'NA'
                ]
            }
            [rec.group, JsonOutput.toJson(samples)]
        }

        ch_chipbinner_inputs = CHIPBINNER_BINS.out.bins.join(ch_chipbinner_samples)
            .map { group, bins_file, samples_json -> [bins_file, samples_json, group] }

        CHIPBINNER_COUNTS(
            ch_chipbinner_inputs.map { bins, samples_json, group -> bins },
            ch_chipbinner_inputs.map { bins, samples_json, group -> samples_json },
            ch_chipbinner_inputs.map { bins, samples_json, group -> group },
            ch_use_spikein,
            Channel.value(params.chipbinner_pseudocount)
        )
        ch_versions = ch_versions.mix(CHIPBINNER_COUNTS.out.versions)

        CHIPBINNER_HDBSCAN_GRID(
            CHIPBINNER_COUNTS.out.normalized,
            ch_chipbinner_inputs.map { bins, samples_json, group -> group },
            Channel.value(params.chipbinner_hdbscan_grid_min_cluster_size),
            Channel.value(params.chipbinner_hdbscan_grid_min_samples)
        )
        ch_versions = ch_versions.mix(CHIPBINNER_HDBSCAN_GRID.out.versions)

        CHIPBINNER_ROTS(
            CHIPBINNER_COUNTS.out.normalized,
            CHIPBINNER_HDBSCAN_GRID.out.clusters,
            ch_chipbinner_inputs.map { bins, samples_json, group -> samples_json },
            Channel.value(treated),
            Channel.value(control),
            ch_chipbinner_inputs.map { bins, samples_json, group -> group },
            Channel.value(params.chipbinner_fdr),
            Channel.value(params.chipbinner_lfc),
            Channel.value(params.chipbinner_bootstrap),
            Channel.value(params.chipbinner_k_value)
        )
        ch_versions = ch_versions.mix(CHIPBINNER_ROTS.out.versions)
        ch_summary_files = ch_summary_files.mix(CHIPBINNER_ROTS.out.summary.map { group, file -> file })
        ch_chipbinner_significant = CHIPBINNER_ROTS.out.up

        if (params.differential_annotate) {
            ch_chipbinner_annotate = CHIPBINNER_ROTS.out.results.map { group, file ->
                def meta = [publish_dir: "${params.outdir}/03_peak_calling/06_differential/02_chipbinner/${group}/differential"]
                [meta, file]
            }
            ANNOTATE_REGIONS(
                ch_chipbinner_annotate,
                ch_gene_bed,
                ch_gtf
            )
            ch_versions = ch_versions.mix(ANNOTATE_REGIONS.out.versions)
        }
    }

    // SPAN differential per group
    if (params.run_span_diff) {
        if (!params.omnipeaks_jar) {
            exit 1, "--omnipeaks_jar is required when --run_span_diff is enabled"
        }
        SPAN_CAPABILITY_PROBE(Channel.value(file(params.omnipeaks_jar)))
        ch_versions = ch_versions.mix(SPAN_CAPABILITY_PROBE.out.versions)

        ch_span_caps = SPAN_CAPABILITY_PROBE.out.capabilities
            .map { file -> new JsonSlurper().parse(file) }
        ch_span_has_compare = ch_span_caps.map { it.has_compare ?: false }

        ch_span_groups = ch_valid_groups.map { rec ->
            def treated_bams = rec.records.findAll { it.condition == treated }.collect { it.bam }
            def control_bams = rec.records.findAll { it.condition == control }.collect { it.bam }
            def meta = [id: rec.group, group: rec.group, treated: treated, control: control]
            [meta, treated_bams, control_bams]
        }

        def span_mode = params.span_diff_mode?.toString()?.toLowerCase() ?: 'auto'

        ch_span_native_groups = Channel.empty()
        if (span_mode == 'native') {
            ch_span_native_groups = ch_span_groups
        } else if (span_mode == 'auto') {
            ch_span_native_groups = ch_span_groups.combine(ch_span_has_compare)
                .filter { meta, treated_bams, control_bams, has_compare -> has_compare }
                .map { meta, treated_bams, control_bams, has_compare -> [meta, treated_bams, control_bams] }
        }

        if (span_mode == 'native' || span_mode == 'auto') {
            SPAN_COMPARE(
                ch_span_native_groups,
            ch_chrom_sizes_single,
            Channel.value(file(params.omnipeaks_jar)),
            Channel.value(params.span_diff_gap),
            Channel.value(params.span_diff_bin),
            Channel.value(params.span_diff_fdr),
            Channel.value(params.span_diff_java_heap)
            )
            ch_versions = ch_versions.mix(SPAN_COMPARE.out.versions)
            ch_summary_files = ch_summary_files.mix(SPAN_COMPARE.out.summary.map { meta, file -> file })
            ch_span_significant = ch_span_significant.mix(SPAN_COMPARE.out.up.map { meta, file -> [meta.group, file] })

            if (params.differential_annotate) {
                ch_span_annotate = SPAN_COMPARE.out.tsv.map { meta, file ->
                    def publish_meta = [publish_dir: "${params.outdir}/03_peak_calling/06_differential/03_span/${meta.group}"]
                    [publish_meta, file]
                }
                ANNOTATE_REGIONS(
                    ch_span_annotate,
                    ch_gene_bed,
                    ch_gtf
                )
                ch_versions = ch_versions.mix(ANNOTATE_REGIONS.out.versions)
            }
        }

        ch_span_fallback_groups = Channel.empty()
        if (span_mode == 'fallback') {
            ch_span_fallback_groups = ch_valid_groups
        } else if (span_mode == 'auto') {
            ch_span_fallback_groups = ch_valid_groups.combine(ch_span_has_compare)
                .filter { rec, has_compare -> !has_compare }
                .map { rec, has_compare -> rec }
        }

        if (span_mode == 'fallback' || span_mode == 'auto') {
            def span_callers = callers.findAll { it.startsWith('span_') }
            if (!span_callers) {
                ch_skipped_records = ch_skipped_records.mix(
                    ch_span_fallback_groups.map { rec -> [method: 'span', group: rec.group, caller: 'NA', reason: 'no_span_peaks', details: 'SPAN caller not in --peakcaller'] }
                )
            } else {
                def span_caller = span_callers[0]
                ch_span_peaks = ch_peaks
                    .filter { meta, peaks -> meta.is_control == false && meta.caller == span_caller }
                    .map { meta, peaks -> [meta.sample_id ?: meta.id, peaks] }
                    .toList()
                    .map { list -> list.collectEntries { [(it[0]): it[1]] } }

                ch_span_samples = ch_span_fallback_groups.combine(ch_span_peaks)
                    .map { rec, peaks_map ->
                        def missing = rec.records.findAll { !peaks_map.containsKey(it.sample_id) }.collect { it.sample_id }
                        def samples = rec.records.collect { row ->
                            [
                                sample_id: row.sample_id,
                                group: row.group,
                                condition: row.condition,
                                replicate: row.replicate,
                                bam: row.bam_path,
                                peaks: peaks_map.get(row.sample_id)?.toString(),
                                caller: span_caller,
                                spikein_scale_factor: row.spikein_scale_factor ?: 'NA'
                            ]
                        }
                        [group: rec.group, samples: samples, missing: missing]
                    }

                ch_span_samples.branch { valid: it.missing.isEmpty(); invalid: !it.missing.isEmpty() }.set { ch_span_split }
                ch_span_valid = ch_span_split.valid
                ch_span_invalid = ch_span_split.invalid

                ch_skipped_records = ch_skipped_records.mix(
                    ch_span_invalid.map { rec -> [method: 'span', group: rec.group, caller: span_caller, reason: 'missing_peaks', details: rec.missing.join(',')] }
                )

                MAKE_SPAN_SAMPLESHEET(
                    ch_span_valid.map { rec -> JsonOutput.toJson(rec.samples) },
                    ch_span_valid.map { rec -> rec.group },
                    Channel.value('span_fallback')
                )
                ch_versions = ch_versions.mix(MAKE_SPAN_SAMPLESHEET.out.versions)

                SPAN_FALLBACK_DIFF(
                    MAKE_SPAN_SAMPLESHEET.out.samplesheet.map { group, caller, sheet -> sheet },
                    Channel.value(treated),
                    Channel.value(control),
                    MAKE_SPAN_SAMPLESHEET.out.samplesheet.map { group, caller, sheet -> group },
                    ch_use_spikein,
                    Channel.value(params.span_diff_fdr),
                    Channel.value(params.diffbind_lfc),
                    Channel.value(params.diffbind_min_overlap),
                    Channel.value(params.diffbind_summits),
                    Channel.value(params.diffbind_backend)
                )
                ch_versions = ch_versions.mix(SPAN_FALLBACK_DIFF.out.versions)
                ch_summary_files = ch_summary_files.mix(SPAN_FALLBACK_DIFF.out.summary.map { group, file -> file })
                ch_span_significant = ch_span_significant.mix(SPAN_FALLBACK_DIFF.out.up.map { group, file -> [group, file] })

                if (params.differential_annotate) {
                    ch_span_fallback_annotate = SPAN_FALLBACK_DIFF.out.results.map { group, file ->
                        def meta = [publish_dir: "${params.outdir}/03_peak_calling/06_differential/03_span/${group}"]
                        [meta, file]
                    }
                    ANNOTATE_REGIONS(
                        ch_span_fallback_annotate,
                        ch_gene_bed,
                        ch_gtf
                    )
                    ch_versions = ch_versions.mix(ANNOTATE_REGIONS.out.versions)
                }
            }
        }
    }

    if (params.differential_cross_compare) {
        def primary_caller = callers ? callers[0] : null
        diffbind_beds = ch_diffbind_significant.map { group, caller, file ->
            [group: group, caller: caller, path: file.toString()]
        }
        method_beds = Channel.empty()
        if (primary_caller) {
            method_beds = method_beds.mix(
                ch_diffbind_significant
                    .filter { group, caller, file -> caller == primary_caller }
                    .map { group, caller, file -> [group: group, method: 'diffbind', caller: caller, path: file.toString()] }
            )
        }
        method_beds = method_beds.mix(ch_chipbinner_significant.map { group, file -> [group: group, method: 'chipbinner', caller: 'NA', path: file.toString()] })
        method_beds = method_beds.mix(ch_span_significant.map { group, file -> [group: group, method: 'span', caller: 'NA', path: file.toString()] })

        callers_tsv = diffbind_beds
            .map { rec -> [rec.group, rec.caller, rec.path].join('\t') }
            .toSortedList()
            .map { lines -> lines.join('\n') + '\n' }
            .collectFile(name: 'overlap_callers.tsv')
        methods_tsv = method_beds
            .map { rec -> [rec.group, rec.method, rec.caller ?: 'NA', rec.path].join('\t') }
            .toSortedList()
            .map { lines -> lines.join('\n') + '\n' }
            .collectFile(name: 'overlap_methods.tsv')

        DIFFERENTIAL_OVERLAP(callers_tsv, methods_tsv)
        ch_versions = ch_versions.mix(DIFFERENTIAL_OVERLAP.out.versions)
    }

    ch_summary_out = Channel.empty()
    ch_skipped_out = Channel.empty()
    if (params.differential_run_multiqc) {
        DIFFERENTIAL_SUMMARY_MERGE(
            ch_summary_files.collect().ifEmpty([]),
            ch_skipped_records.toList().map { list -> JsonOutput.toJson(list) },
            Channel.value(file("${projectDir}/assets/multiqc/differential_summary_header.txt")),
            Channel.value(file("${projectDir}/assets/multiqc/differential_skipped_header.txt"))
        )
        ch_versions = ch_versions.mix(DIFFERENTIAL_SUMMARY_MERGE.out.versions)
        ch_summary_out = DIFFERENTIAL_SUMMARY_MERGE.out.summary
        ch_skipped_out = DIFFERENTIAL_SUMMARY_MERGE.out.skipped
    }

    emit:
    summary = ch_summary_out
    skipped = ch_skipped_out
    versions = ch_versions
}
