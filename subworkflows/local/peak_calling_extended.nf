/*
 * Peak calling with condition-aware controls and extended caller support
 */

include { SEACR_CALLPEAK } from "../../modules/nf-core/seacr/callpeak/main"
include { MACS2_CALLPEAK as MACS2_CALLPEAK_LEGACY } from "../../modules/nf-core/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_CALLPEAK_NARROW } from "../../modules/nf-core/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_CALLPEAK_BROAD } from "../../modules/nf-core/macs2/callpeak/main"
include { GOPEAKS_CALLPEAK as GOPEAKS_CALLPEAK_NARROW } from "../../modules/local/gopeaks_callpeak"
include { GOPEAKS_CALLPEAK as GOPEAKS_CALLPEAK_BROAD } from "../../modules/local/gopeaks_callpeak"
include { EPIC2_CALLPEAK as EPIC2_CALLPEAK_200 } from "../../modules/local/epic2_callpeak"
include { EPIC2_CALLPEAK as EPIC2_CALLPEAK_150 } from "../../modules/local/epic2_callpeak"
include { EPIC2_CALLPEAK as EPIC2_CALLPEAK_25 } from "../../modules/local/epic2_callpeak"
include { SPAN_OMNIPEAKS_ANALYZE as SPAN_OMNIPEAKS_DEFAULT } from "../../modules/local/span_omnipeaks_analyze"
include { SPAN_OMNIPEAKS_ANALYZE as SPAN_OMNIPEAKS_STRINGENT } from "../../modules/local/span_omnipeaks_analyze"
include { POOL_IGG_CONTROLS } from "../../modules/local/pool_igg_controls"

workflow PEAK_CALLING_EXTENDED {
    take:
    bedgraph_target   // channel: [ val(meta), bedgraph ]
    bedgraph_control  // channel: [ val(meta), bedgraph ]
    bam_target        // channel: [ val(meta), bam ]
    bam_control       // channel: [ val(meta), bam ]
    chrom_sizes       // channel: [ chrom_sizes ]
    callers           // list of caller ids

    main:
    ch_versions = Channel.empty()
    ch_peaks_all = Channel.empty()
    ch_macs2_summits = Channel.empty()
    ch_control_fallbacks = Channel.empty()
    ch_pooled_controls_bam = Channel.empty()
    ch_pooled_controls_bai = Channel.empty()
    ch_gopeaks_json = Channel.empty()

    def primary_caller = callers ? callers[0] : 'seacr'
    def control_required_callers = callers.findAll { it.startsWith('epic2_') || it.startsWith('span_') }
    if (control_required_callers && !params.use_control) {
        exit 1, "Peak callers requiring controls (${control_required_callers.join(', ')}) cannot be run with --use_control false. Remove these callers or enable controls with --use_control true."
    }

    // Build control condition map by control group
    def build_control_condition_map = { ch_control ->
        ch_control
            .map { meta, file -> [meta.group, meta.condition] }
            .toList()
            .map { list ->
                def control_map = [:].withDefault { [] }
                list.each { entry ->
                    def group = entry[0]
                    def condition = entry[1]
                    control_map[group] = (control_map[group] ?: []) + [condition]
                }
                control_map
            }
    }

    // Annotate target metas with control_condition (preferred condition, fallback to NA or first available)
    def add_control_condition = { ch_target, ch_control ->
        def control_map_ch = build_control_condition_map(ch_control)
        ch_target
            .combine(control_map_ch)
            .map { meta, file, control_map ->
                def conditions = control_map.get(meta.control_group, [])
                def control_condition = meta.condition
                if (conditions) {
                    if (conditions.contains(meta.condition)) {
                        control_condition = meta.condition
                    } else if (conditions.contains('NA')) {
                        control_condition = 'NA'
                    } else {
                        control_condition = conditions.sort()[0]
                    }
                }
                def meta_out = meta + [control_condition: control_condition]
                [meta_out, file]
            }
    }

    // Pair targets with controls by control_group + control_condition and replicate logic
    // Never drop targets when controls are missing (return null control instead)
    def pair_targets_controls = { ch_target, ch_control ->
        def control_map_ch = ch_control
            .map { meta, file -> ["${meta.group}__${meta.condition}", [meta, file]] }
            .groupTuple(by: [0])
            .toList()
            .map { list ->
                def map = [:]
                (list ?: []).each { entry ->
                    map[entry[0]] = entry[1]
                }
                map
            }

        def target_count_ch = ch_target
            .map { meta, file -> ["${meta.control_group}__${meta.control_condition}", meta.replicate] }
            .groupTuple(by: [0])
            .map { key, reps -> [key, reps.size()] }
            .toList()
            .map { list ->
                def map = [:]
                (list ?: []).each { entry ->
                    map[entry[0]] = entry[1]
                }
                map
            }

        ch_target
            .combine(control_map_ch)
            .combine(target_count_ch)
            .flatMap { meta, file, control_map, target_count_map ->
                def key = "${meta.control_group}__${meta.control_condition}"
                def control_list = control_map.get(key)
                if (!control_list) {
                    return [[meta + [control_missing: true], file, null]]
                }
                def control_by_rep = control_list.collectEntries { [(it[0].replicate): it[1]] }
                def control_default = control_list.sort { it[0].replicate }[0][1]
                def control_reps_count = control_list.size()
                def target_reps_count = target_count_map.get(key) ?: 1
                def control_file = (control_reps_count == target_reps_count && control_by_rep.containsKey(meta.replicate)) ?
                    control_by_rep[meta.replicate] : control_default
                return [[meta + [control_missing: false], file, control_file]]
            }
    }

    // Annotate targets with control condition
    ch_bedgraph_target_cc = add_control_condition(bedgraph_target, bedgraph_control)
    ch_bam_target_cc = add_control_condition(bam_target, bam_control)

    /*
     * SEACR
     */
    if ('seacr' in callers) {
        if (params.use_control) {
            ch_seacr_pairs = pair_targets_controls(
                ch_bedgraph_target_cc.map { meta, bed -> [meta + [caller: 'seacr'], bed] },
                bedgraph_control
            )
            def ch_seacr_inputs = ch_seacr_pairs.map { meta, bed, control ->
                if (!control) {
                    log.warn "No control found for group '${meta.control_group}' (condition '${meta.control_condition}') - running SEACR without control for ${meta.sample_id ?: meta.id}"
                    return [meta, bed, []]
                }
                [meta, bed, control]
            }
            SEACR_CALLPEAK (
                ch_seacr_inputs,
                params.seacr_peak_threshold
            )
            ch_peaks_all = ch_peaks_all.mix(SEACR_CALLPEAK.out.bed)
            ch_versions = ch_versions.mix(SEACR_CALLPEAK.out.versions)
        } else {
            ch_seacr_nocontrol = ch_bedgraph_target_cc
                .map { meta, bed -> [meta + [caller: 'seacr'], bed, []] }
            SEACR_CALLPEAK (
                ch_seacr_nocontrol,
                params.seacr_peak_threshold
            )
            ch_peaks_all = ch_peaks_all.mix(SEACR_CALLPEAK.out.bed)
            ch_versions = ch_versions.mix(SEACR_CALLPEAK.out.versions)
        }
    }

    /*
     * Legacy MACS2 (with optional control)
     */
    if ('macs2' in callers) {
        if (params.use_control) {
            ch_macs_pairs = pair_targets_controls(
                ch_bam_target_cc.map { meta, bam -> [meta + [caller: 'macs2'], bam] },
                bam_control
            )
            def ch_macs_inputs = ch_macs_pairs.map { meta, bam, control ->
                if (!control) {
                    log.warn "No control found for group '${meta.control_group}' (condition '${meta.control_condition}') - running MACS2 without control for ${meta.sample_id ?: meta.id}"
                    return [meta, bam, []]
                }
                [meta, bam, control]
            }
            MACS2_CALLPEAK_LEGACY (
                ch_macs_inputs,
                params.macs_gsize
            )
            ch_peaks_all = ch_peaks_all.mix(MACS2_CALLPEAK_LEGACY.out.peak)
            ch_macs2_summits = ch_macs2_summits.mix(MACS2_CALLPEAK_LEGACY.out.bed)
            ch_versions = ch_versions.mix(MACS2_CALLPEAK_LEGACY.out.versions)
        } else {
            ch_macs_nocontrol = ch_bam_target_cc
                .map { meta, bam -> [meta + [caller: 'macs2'], bam, []] }
            MACS2_CALLPEAK_LEGACY (
                ch_macs_nocontrol,
                params.macs_gsize
            )
            ch_peaks_all = ch_peaks_all.mix(MACS2_CALLPEAK_LEGACY.out.peak)
            ch_macs2_summits = ch_macs2_summits.mix(MACS2_CALLPEAK_LEGACY.out.bed)
            ch_versions = ch_versions.mix(MACS2_CALLPEAK_LEGACY.out.versions)
        }
    }

    /*
     * MACS2 narrow + broad (no control)
     */
    if ('macs2_narrow' in callers) {
        ch_macs2_narrow = ch_bam_target_cc.map { meta, bam -> [meta + [caller: 'macs2_narrow'], bam, []] }
        MACS2_CALLPEAK_NARROW (
            ch_macs2_narrow,
            params.macs_gsize
        )
        ch_peaks_all = ch_peaks_all.mix(MACS2_CALLPEAK_NARROW.out.peak)
        ch_macs2_summits = ch_macs2_summits.mix(MACS2_CALLPEAK_NARROW.out.bed)
        ch_versions = ch_versions.mix(MACS2_CALLPEAK_NARROW.out.versions)
    }

    if ('macs2_broad' in callers) {
        ch_macs2_broad = ch_bam_target_cc.map { meta, bam -> [meta + [caller: 'macs2_broad'], bam, []] }
        MACS2_CALLPEAK_BROAD (
            ch_macs2_broad,
            params.macs_gsize
        )
        ch_peaks_all = ch_peaks_all.mix(MACS2_CALLPEAK_BROAD.out.peak)
        ch_macs2_summits = ch_macs2_summits.mix(MACS2_CALLPEAK_BROAD.out.bed)
        ch_versions = ch_versions.mix(MACS2_CALLPEAK_BROAD.out.versions)
    }

    /*
     * GoPeaks
     */
    if ('gopeaks_narrow' in callers) {
        GOPEAKS_CALLPEAK_NARROW (
            ch_bam_target_cc.map { meta, bam -> [meta + [caller: 'gopeaks_narrow'], bam] },
            false
        )
        ch_peaks_all = ch_peaks_all.mix(GOPEAKS_CALLPEAK_NARROW.out.peaks)
        ch_gopeaks_json = ch_gopeaks_json.mix(GOPEAKS_CALLPEAK_NARROW.out.json)
        ch_versions = ch_versions.mix(GOPEAKS_CALLPEAK_NARROW.out.versions)
    }

    if ('gopeaks_broad' in callers) {
        GOPEAKS_CALLPEAK_BROAD (
            ch_bam_target_cc.map { meta, bam -> [meta + [caller: 'gopeaks_broad'], bam] },
            true
        )
        ch_peaks_all = ch_peaks_all.mix(GOPEAKS_CALLPEAK_BROAD.out.peaks)
        ch_gopeaks_json = ch_gopeaks_json.mix(GOPEAKS_CALLPEAK_BROAD.out.json)
        ch_versions = ch_versions.mix(GOPEAKS_CALLPEAK_BROAD.out.versions)
    }

    /*
     * Pooled controls for epic2/SPAN
     */
    def needs_pooled_controls = params.use_control && callers.any { it.startsWith('epic2_') || it.startsWith('span_') }
    if (needs_pooled_controls) {
    ch_pooled_inputs = bam_control
        .map { meta, bam -> ["${meta.group}__${meta.condition}", [meta, bam]] }
        .groupTuple(by: [0])
        .map { key, entries ->
                def meta0 = entries[0][0]
                def pooled_meta = meta0 + [
                    id: "${meta0.group}_${meta0.condition}",
                    control_group: meta0.group,
                    is_control: true
                ]
                [pooled_meta, entries.collect { it[1] }]
            }

        POOL_IGG_CONTROLS ( ch_pooled_inputs )
        ch_pooled_controls_bam = POOL_IGG_CONTROLS.out.bam
        ch_pooled_controls_bai = POOL_IGG_CONTROLS.out.bai
        ch_versions = ch_versions.mix(POOL_IGG_CONTROLS.out.versions)
    }

    ch_pooled_pairs = Channel.empty()
    if (needs_pooled_controls) {
        pooled_map_ch = ch_pooled_controls_bam
            .map { meta, bam -> [meta.control_group ?: meta.group, meta.condition, bam] }
            .toList()
            .map { list ->
                def map = [:].withDefault { [] }
                list.each { entry ->
                    def group = entry[0]
                    def condition = entry[1]
                    def bam = entry[2]
                    map[group] = (map[group] ?: []) + [[condition, bam]]
                }
                map
            }

        ch_pooled_pairs = ch_bam_target_cc
            .combine(pooled_map_ch)
            .map { meta, bam, pooled_map ->
                def entries = pooled_map.get(meta.control_group, [])
                if (!entries) {
                    return [meta, bam, null, '', 'missing_control', 'skipped', 'no_control_for_group', null]
                }
                def desired_condition = meta.control_condition ?: meta.condition
                def selected_entry = entries.find { it[0] == desired_condition }
                def used_condition = null
                def control_bam = null
                def status = null
                def reason = null
                if (selected_entry) {
                    used_condition = selected_entry[0]
                    control_bam = selected_entry[1]
                } else {
                    def na_entry = entries.find { it[0] == 'NA' }
                    if (na_entry) {
                        used_condition = na_entry[0]
                        control_bam = na_entry[1]
                    } else {
                        def fallback = entries.sort { it[0] }[0]
                        used_condition = fallback[0]
                        control_bam = fallback[1]
                    }
                }
                def exact_match = (used_condition == meta.condition)
                if (exact_match) {
                    status = 'exact_match'
                    reason = 'exact_condition'
                } else if (used_condition == 'NA') {
                    status = 'fallback_other_condition'
                    reason = 'legacy_na_control'
                } else {
                    status = 'fallback_other_condition'
                    reason = 'condition_fallback'
                }
                return [meta, bam, control_bam, used_condition, status, 'used', reason, control_bam]
            }
    }

    chrom_sizes_val = chrom_sizes
        .collect()
        .map { it instanceof List ? it[0] : it }

    omnipeaks_jar_val = Channel.empty()
    if (params.omnipeaks_jar) {
        omnipeaks_jar_val = Channel.fromPath(params.omnipeaks_jar)
            .collect()
            .map { it instanceof List ? it[0] : it }
    }

    /*
     * epic2
     */
    if ('epic2_200bp' in callers) {
        ch_epic2_200 = ch_pooled_pairs
            .map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                [meta + [caller: 'epic2_200bp'], bam, control_bam, used_condition, status, action, reason, pooled_path]
            }
        ch_epic2_200_branch = ch_epic2_200.branch {
            run: it[2] != null
            skip: it[2] == null
        }
        EPIC2_CALLPEAK_200 (
            ch_epic2_200_branch.run.map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                [meta, bam, control_bam]
            },
            params.epic2_genome,
            200,
            3,
            params.epic2_fdr
        )
        ch_peaks_all = ch_peaks_all.mix(EPIC2_CALLPEAK_200.out.peaks)
        ch_versions = ch_versions.mix(EPIC2_CALLPEAK_200.out.versions)
        ch_control_fallbacks = ch_control_fallbacks.mix(
            ch_epic2_200_branch.run
                .filter { meta, bam, control_bam, used_condition, status, action, reason, pooled_path -> status != 'exact_match' }
                .map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                    [
                        sample_id: meta.sample_id ?: meta.id,
                        group: meta.group,
                        condition: meta.condition,
                        caller_id: 'epic2_200bp',
                        control_group: meta.control_group,
                        selected_control_condition: used_condition,
                        status: status,
                        action: action,
                        reason: reason,
                        pooled_control_path: pooled_path
                    ]
                }
        )
        ch_control_fallbacks = ch_control_fallbacks.mix(
            ch_epic2_200_branch.skip
                .map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                    log.warn("Skipping epic2_200bp for sample ${meta.sample_id ?: meta.id} (control_group=${meta.control_group}, condition=${meta.condition}) - no pooled control available.")
                    [
                        sample_id: meta.sample_id ?: meta.id,
                        group: meta.group,
                        condition: meta.condition,
                        caller_id: 'epic2_200bp',
                        control_group: meta.control_group,
                        selected_control_condition: used_condition ?: '',
                        status: status ?: 'missing_control',
                        action: 'skipped',
                        reason: reason ?: 'missing_control',
                        pooled_control_path: pooled_path ?: ''
                    ]
                }
        )
    }

    if ('epic2_150bp' in callers) {
        ch_epic2_150 = ch_pooled_pairs
            .map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                [meta + [caller: 'epic2_150bp'], bam, control_bam, used_condition, status, action, reason, pooled_path]
            }
        ch_epic2_150_branch = ch_epic2_150.branch {
            run: it[2] != null
            skip: it[2] == null
        }
        EPIC2_CALLPEAK_150 (
            ch_epic2_150_branch.run.map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                [meta, bam, control_bam]
            },
            params.epic2_genome,
            150,
            2,
            params.epic2_fdr
        )
        ch_peaks_all = ch_peaks_all.mix(EPIC2_CALLPEAK_150.out.peaks)
        ch_versions = ch_versions.mix(EPIC2_CALLPEAK_150.out.versions)
        ch_control_fallbacks = ch_control_fallbacks.mix(
            ch_epic2_150_branch.run
                .filter { meta, bam, control_bam, used_condition, status, action, reason, pooled_path -> status != 'exact_match' }
                .map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                    [
                        sample_id: meta.sample_id ?: meta.id,
                        group: meta.group,
                        condition: meta.condition,
                        caller_id: 'epic2_150bp',
                        control_group: meta.control_group,
                        selected_control_condition: used_condition,
                        status: status,
                        action: action,
                        reason: reason,
                        pooled_control_path: pooled_path
                    ]
                }
        )
        ch_control_fallbacks = ch_control_fallbacks.mix(
            ch_epic2_150_branch.skip
                .map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                    log.warn("Skipping epic2_150bp for sample ${meta.sample_id ?: meta.id} (control_group=${meta.control_group}, condition=${meta.condition}) - no pooled control available.")
                    [
                        sample_id: meta.sample_id ?: meta.id,
                        group: meta.group,
                        condition: meta.condition,
                        caller_id: 'epic2_150bp',
                        control_group: meta.control_group,
                        selected_control_condition: used_condition ?: '',
                        status: status ?: 'missing_control',
                        action: 'skipped',
                        reason: reason ?: 'missing_control',
                        pooled_control_path: pooled_path ?: ''
                    ]
                }
        )
    }

    if ('epic2_25bp' in callers) {
        ch_epic2_25 = ch_pooled_pairs
            .map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                [meta + [caller: 'epic2_25bp'], bam, control_bam, used_condition, status, action, reason, pooled_path]
            }
        ch_epic2_25_branch = ch_epic2_25.branch {
            run: it[2] != null
            skip: it[2] == null
        }
        EPIC2_CALLPEAK_25 (
            ch_epic2_25_branch.run.map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                [meta, bam, control_bam]
            },
            params.epic2_genome,
            25,
            2,
            params.epic2_fdr
        )
        ch_peaks_all = ch_peaks_all.mix(EPIC2_CALLPEAK_25.out.peaks)
        ch_versions = ch_versions.mix(EPIC2_CALLPEAK_25.out.versions)
        ch_control_fallbacks = ch_control_fallbacks.mix(
            ch_epic2_25_branch.run
                .filter { meta, bam, control_bam, used_condition, status, action, reason, pooled_path -> status != 'exact_match' }
                .map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                    [
                        sample_id: meta.sample_id ?: meta.id,
                        group: meta.group,
                        condition: meta.condition,
                        caller_id: 'epic2_25bp',
                        control_group: meta.control_group,
                        selected_control_condition: used_condition,
                        status: status,
                        action: action,
                        reason: reason,
                        pooled_control_path: pooled_path
                    ]
                }
        )
        ch_control_fallbacks = ch_control_fallbacks.mix(
            ch_epic2_25_branch.skip
                .map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                    log.warn("Skipping epic2_25bp for sample ${meta.sample_id ?: meta.id} (control_group=${meta.control_group}, condition=${meta.condition}) - no pooled control available.")
                    [
                        sample_id: meta.sample_id ?: meta.id,
                        group: meta.group,
                        condition: meta.condition,
                        caller_id: 'epic2_25bp',
                        control_group: meta.control_group,
                        selected_control_condition: used_condition ?: '',
                        status: status ?: 'missing_control',
                        action: 'skipped',
                        reason: reason ?: 'missing_control',
                        pooled_control_path: pooled_path ?: ''
                    ]
                }
        )
    }

    /*
     * SPAN / OMNIPEAKS
     */
    if ('span_default' in callers) {
        ch_span_default = ch_pooled_pairs
            .map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                [meta + [caller: 'span_default'], bam, control_bam, used_condition, status, action, reason, pooled_path]
            }
        ch_span_default_branch = ch_span_default.branch {
            run: it[2] != null
            skip: it[2] == null
        }
        SPAN_OMNIPEAKS_DEFAULT (
            ch_span_default_branch.run.map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                [meta, bam, control_bam]
            },
            chrom_sizes_val,
            omnipeaks_jar_val,
            5,
            '',
            params.span_java_heap
        )
        ch_peaks_all = ch_peaks_all.mix(SPAN_OMNIPEAKS_DEFAULT.out.peaks)
        ch_versions = ch_versions.mix(SPAN_OMNIPEAKS_DEFAULT.out.versions)
        ch_control_fallbacks = ch_control_fallbacks.mix(
            ch_span_default_branch.run
                .filter { meta, bam, control_bam, used_condition, status, action, reason, pooled_path -> status != 'exact_match' }
                .map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                    [
                        sample_id: meta.sample_id ?: meta.id,
                        group: meta.group,
                        condition: meta.condition,
                        caller_id: 'span_default',
                        control_group: meta.control_group,
                        selected_control_condition: used_condition,
                        status: status,
                        action: action,
                        reason: reason,
                        pooled_control_path: pooled_path
                    ]
                }
        )
        ch_control_fallbacks = ch_control_fallbacks.mix(
            ch_span_default_branch.skip
                .map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                    log.warn("Skipping span_default for sample ${meta.sample_id ?: meta.id} (control_group=${meta.control_group}, condition=${meta.condition}) - no pooled control available.")
                    [
                        sample_id: meta.sample_id ?: meta.id,
                        group: meta.group,
                        condition: meta.condition,
                        caller_id: 'span_default',
                        control_group: meta.control_group,
                        selected_control_condition: used_condition ?: '',
                        status: status ?: 'missing_control',
                        action: 'skipped',
                        reason: reason ?: 'missing_control',
                        pooled_control_path: pooled_path ?: ''
                    ]
                }
        )
    }

    if ('span_stringent' in callers) {
        ch_span_stringent = ch_pooled_pairs
            .map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                [meta + [caller: 'span_stringent'], bam, control_bam, used_condition, status, action, reason, pooled_path]
            }
        ch_span_stringent_branch = ch_span_stringent.branch {
            run: it[2] != null
            skip: it[2] == null
        }
        SPAN_OMNIPEAKS_STRINGENT (
            ch_span_stringent_branch.run.map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                [meta, bam, control_bam]
            },
            chrom_sizes_val,
            omnipeaks_jar_val,
            2,
            params.span_stringent_fdr,
            params.span_java_heap
        )
        ch_peaks_all = ch_peaks_all.mix(SPAN_OMNIPEAKS_STRINGENT.out.peaks)
        ch_versions = ch_versions.mix(SPAN_OMNIPEAKS_STRINGENT.out.versions)
        ch_control_fallbacks = ch_control_fallbacks.mix(
            ch_span_stringent_branch.run
                .filter { meta, bam, control_bam, used_condition, status, action, reason, pooled_path -> status != 'exact_match' }
                .map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                    [
                        sample_id: meta.sample_id ?: meta.id,
                        group: meta.group,
                        condition: meta.condition,
                        caller_id: 'span_stringent',
                        control_group: meta.control_group,
                        selected_control_condition: used_condition,
                        status: status,
                        action: action,
                        reason: reason,
                        pooled_control_path: pooled_path
                    ]
                }
        )
        ch_control_fallbacks = ch_control_fallbacks.mix(
            ch_span_stringent_branch.skip
                .map { meta, bam, control_bam, used_condition, status, action, reason, pooled_path ->
                    log.warn("Skipping span_stringent for sample ${meta.sample_id ?: meta.id} (control_group=${meta.control_group}, condition=${meta.condition}) - no pooled control available.")
                    [
                        sample_id: meta.sample_id ?: meta.id,
                        group: meta.group,
                        condition: meta.condition,
                        caller_id: 'span_stringent',
                        control_group: meta.control_group,
                        selected_control_condition: used_condition ?: '',
                        status: status ?: 'missing_control',
                        action: 'skipped',
                        reason: reason ?: 'missing_control',
                        pooled_control_path: pooled_path ?: ''
                    ]
                }
        )
    }

    // Identify primary and secondary peaks by caller id
    ch_peaks_all
        .branch {
            primary: it[0].caller == primary_caller
            secondary: it[0].caller != primary_caller
        }
        .set { ch_peaks_split }

    ch_peaks_primary = ch_peaks_split.primary
    ch_peaks_secondary = ch_peaks_split.secondary

    emit:
    peaks = ch_peaks_all
    peaks_primary = ch_peaks_primary
    peaks_secondary = ch_peaks_secondary
    macs2_summits = ch_macs2_summits
    pooled_controls_bam = ch_pooled_controls_bam
    pooled_controls_bai = ch_pooled_controls_bai
    control_fallbacks = ch_control_fallbacks
    gopeaks_json = ch_gopeaks_json
    versions = ch_versions
}
