/*
 * Calculate Peak-based metrics and QC
*/

include { PEAK_FRIP                            } from "../../modules/local/peak_frip"
include { PEAK_COUNTS as PRIMARY_PEAK_COUNTS   } from "../../modules/local/peak_counts"
include { PEAK_COUNTS as CONSENSUS_PEAK_COUNTS } from "../../modules/local/peak_counts"
include { CUT as CUT_CALC_REPROD               } from "../../modules/local/linux/cut"
include { BEDTOOLS_INTERSECT                   } from "../../modules/nf-core/bedtools/intersect/main.nf"
include { CALCULATE_PEAK_REPROD                } from "../../modules/local/python/peak_reprod"
include { PLOT_CONSENSUS_PEAKS                 } from '../../modules/local/python/plot_consensus_peaks'
include { PEAK_QC_TABLE_REPORT as PEAK_QC_FRIP_REPORT      } from "../../modules/local/peak_qc_table_report"
include { PEAK_QC_TABLE_REPORT as PEAK_QC_COUNTS_REPORT    } from "../../modules/local/peak_qc_table_report"
include { PEAK_QC_TABLE_REPORT as PEAK_QC_CONSENSUS_REPORT } from "../../modules/local/peak_qc_table_report"
include { PEAK_QC_TABLE_REPORT as PEAK_QC_REPROD_REPORT    } from "../../modules/local/peak_qc_table_report"

workflow PEAK_QC {
    take:
    peaks                               // channel: [ val(meta), [ bed ] ]
    peaks_with_ids                      // channel: [ val(meta), [ bed ] ]
    consensus_peaks                     // channel: [ val(meta), [ bed ] ]
    consensus_peaks_unfiltered          // channel: [ val(meta), [ bed ] ]
    fragments_bed                       // channel: [ val(meta), [ bed ] ]
    flagstat                            // channel: [ val(meta), [ flagstat ] ]
    min_frip_overlap                    // val
    consensus_grouping                  // val: group or group_condition
    frip_score_header_multiqc           // file
    peak_count_header_multiqc           // file
    peak_count_consensus_header_multiqc // file
    peak_reprod_header_multiqc          // file

    main:
    ch_versions = Channel.empty()

    /*
    * CHANNEL: Combine channel together for frip calculation
    */
    peaks
    .map { row -> [row[0].id, row ].flatten()}
    .join ( fragments_bed.map { row -> [row[0].id, row ].flatten()} )
    .join ( flagstat.map { row -> [row[0].id, row ].flatten()} )
    .map { row -> [ row[1], row[2], row[4], row[6] ]}
    .set { ch_frip }

    /*
    * MODULE: Calculate frip scores for sample peaks
    */
    PEAK_FRIP(
        ch_frip,
        frip_score_header_multiqc,
        min_frip_overlap
    )
    ch_versions = ch_versions.mix(PEAK_FRIP.out.versions)
    // PEAK_FRIP.out.frip_mqc | view

    /*
    * MODULE: Calculate peak counts for sample peaks
    */
    PRIMARY_PEAK_COUNTS(
        peaks,
        peak_count_header_multiqc
    )
    ch_versions = ch_versions.mix(PRIMARY_PEAK_COUNTS.out.versions)
    // PRIMARY_PEAK_COUNTS.out.count_mqc | view

    /*
    * MODULE: Calculate peak counts for consensus peaks
    */
    CONSENSUS_PEAK_COUNTS(
        consensus_peaks,
        peak_count_consensus_header_multiqc
    )
    ch_versions = ch_versions.mix(CONSENSUS_PEAK_COUNTS.out.versions)
    // CONSENSUS_PEAK_COUNTS.out.count_mqc | view

    /*
    * MODULE: Trim unwanted columns for downstream reporting
    */
    CUT_CALC_REPROD (
        peaks_with_ids
    )
    ch_versions = ch_versions.mix(CUT_CALC_REPROD.out.versions)

    /*
    * CHANNEL: Group samples based on group and filter for groups that have more than one file
    */
    CUT_CALC_REPROD.out.file
    .map { row ->
        def group_key = consensus_grouping == 'group_condition' ? row[0].group_condition : row[0].group
        [ "${group_key}__${row[0].caller}", row[1], row[0] ]
    }
    .groupTuple(by: [0])
    .map { row ->
        def meta0 = row[2][0]
        def conditions = row[2].collect { it.condition }.unique()
        def condition_label = consensus_grouping == 'group_condition' ? meta0.condition : (conditions.size() == 1 ? conditions[0] : 'all')
        def group_key = consensus_grouping == 'group_condition' ? meta0.group_condition : meta0.group
        def meta = [id: "${group_key}_${meta0.caller}", group: meta0.group, condition: condition_label, caller: meta0.caller]
        [ meta, row[1].flatten() ]
    }
    .map { row -> [ row[0], row[1], row[1].size() ] }
    .filter { row -> row[2] > 1 }
    .map { row -> [ row[0], row[1] ] }
    .set { ch_peak_bed_group }
    //ch_peak_bed_group | view

    /*
    * CHANNEL: Per group, create a channel per one against all combination
    */
    ch_peak_bed_group.flatMap{
        row ->
        def new_output = []
        row[1].each{ file ->
            def files_copy = row[1].collect()
            files_copy.remove(files_copy.indexOf(file))
            new_output.add([[id: file.name.split("\\.")[0]], file, files_copy])
        }
        new_output
    }
    .set { ch_beds_intersect }
    //EXAMPLE CHANNEL STRUCT: [[META], BED (-a), [BED...n] (-b)]
    //ch_beds_intersect | view

    /*
    * MODULE: Find intra-group overlap
    */
    BEDTOOLS_INTERSECT (
        ch_beds_intersect,
        [[:],[]]
    )
    ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions)
    //EXAMPLE CHANNEL STRUCT: [[META], BED]
    //BEDTOOLS_INTERSECT.out.intersect | view

    /*
    * MODULE: Use overlap to calculate a peak repro %
    */
    CALCULATE_PEAK_REPROD (
        BEDTOOLS_INTERSECT.out.intersect,
        peak_reprod_header_multiqc
    )
    ch_versions = ch_versions.mix(CALCULATE_PEAK_REPROD.out.versions)
    //EXAMPLE CHANNEL STRUCT: [[META], TSV]
    //CALCULATE_PEAK_REPROD.out.tsv

    /*
    * MODULE: Write caller-aware peak QC tables
    */
    def sample_header_prefix = "sample_id\tgroup\tcondition\treplicate\tcaller_id"

    PEAK_FRIP.out.frip_score
    .map { meta, score_file ->
        def score = score_file.text.trim()
        def replicate = meta.replicate ?: 'NA'
        def group = meta.group ?: 'NA'
        def condition = meta.condition ?: 'NA'
        def caller = meta.caller ?: 'NA'
        def sample_id = meta.sample_id ?: meta.id
        [sample_id, group, condition, replicate, caller, score].join('\t')
    }
    .toList()
    .ifEmpty([])
    .map { rows -> ['peak_frip_scores', "${sample_header_prefix}\tfrip_score", rows] }
    .set { ch_frip_table }

    PEAK_QC_FRIP_REPORT(ch_frip_table)
    ch_versions = ch_versions.mix(PEAK_QC_FRIP_REPORT.out.versions)

    PRIMARY_PEAK_COUNTS.out.count_value
    .map { meta, count_file ->
        def count = count_file.text.trim()
        def replicate = meta.replicate ?: 'NA'
        def group = meta.group ?: 'NA'
        def condition = meta.condition ?: 'NA'
        def caller = meta.caller ?: 'NA'
        def sample_id = meta.sample_id ?: meta.id
        [sample_id, group, condition, replicate, caller, count].join('\t')
    }
    .toList()
    .ifEmpty([])
    .map { rows -> ['peak_counts', "${sample_header_prefix}\tpeak_count", rows] }
    .set { ch_peak_count_table }

    PEAK_QC_COUNTS_REPORT(ch_peak_count_table)
    ch_versions = ch_versions.mix(PEAK_QC_COUNTS_REPORT.out.versions)

    CONSENSUS_PEAK_COUNTS.out.count_value
    .map { meta, count_file ->
        def count = count_file.text.trim()
        def group = meta.group ?: 'NA'
        def condition = meta.condition ?: 'NA'
        def caller = meta.caller ?: 'NA'
        def sample_id = meta.sample_id ?: meta.id
        [sample_id, group, condition, 'NA', caller, count].join('\t')
    }
    .toList()
    .ifEmpty([])
    .map { rows -> ['consensus_peak_counts', "${sample_header_prefix}\tconsensus_peak_count", rows] }
    .set { ch_consensus_count_table }

    PEAK_QC_CONSENSUS_REPORT(ch_consensus_count_table)
    ch_versions = ch_versions.mix(PEAK_QC_CONSENSUS_REPORT.out.versions)

    CALCULATE_PEAK_REPROD.out.tsv
    .map { meta, repro_file ->
        def parts = repro_file.text.trim().split('\t')
        def value = parts.size() > 1 ? parts[1] : ''
        def group = meta.group ?: 'NA'
        def condition = meta.condition ?: 'NA'
        def caller = meta.caller ?: 'NA'
        def sample_id = meta.sample_id ?: meta.id
        [sample_id, group, condition, 'NA', caller, value].join('\t')
    }
    .toList()
    .ifEmpty([])
    .map { rows -> ['peak_reproducibility', "${sample_header_prefix}\tpeak_reproducibility_percent", rows] }
    .set { ch_reprod_table }

    PEAK_QC_REPROD_REPORT(ch_reprod_table)
    ch_versions = ch_versions.mix(PEAK_QC_REPROD_REPORT.out.versions)

    /*
    * CHANNEL: Prep for upset input
    */
    consensus_peaks_unfiltered
    .map { row -> [ row[0].caller, row[1] ] }
    .groupTuple(by: [0])
    .map { row ->
        def output = []
        row[1].each{ v -> output.add(v) }
        output
    }
    .set { ch_merged_bed_sorted }

    /*
    * MODULE: Plot upset plots for sample peaks
    */
    PLOT_CONSENSUS_PEAKS (
        ch_merged_bed_sorted.ifEmpty([])
    )
    ch_versions = ch_versions.mix(PLOT_CONSENSUS_PEAKS.out.versions)

    emit:
    primary_frip_mqc    = PEAK_FRIP.out.frip_mqc              // channel: [ val(meta), [ mqc ] ]
    primary_count_mqc   = PRIMARY_PEAK_COUNTS.out.count_mqc   // channel: [ val(meta), [ mqc ] ]
    consensus_count_mqc = CONSENSUS_PEAK_COUNTS.out.count_mqc // channel: [ val(meta), [ mqc ] ]
    reprod_perc_mqc     = CALCULATE_PEAK_REPROD.out.mqc       // channel: [ val(meta), [ mqc ] ]

    versions = ch_versions // channel: [ versions.yml ]
}
