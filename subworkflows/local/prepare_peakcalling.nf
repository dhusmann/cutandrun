/*
 * Convert bam files to bedgraph and bigwig with apropriate normalisation
 */

include { BEDTOOLS_GENOMECOV    } from "../../modules/nf-core/bedtools/genomecov/main"
include { DEEPTOOLS_BAMCOVERAGE } from "../../modules/local/for_patch/deeptools/bamcoverage/main"
include { BEDTOOLS_SORT         } from "../../modules/local/for_patch/bedtools/sort/main"
include { UCSC_BEDCLIP          } from "../../modules/nf-core/ucsc/bedclip/main"
include { UCSC_BEDGRAPHTOBIGWIG } from "../../modules/nf-core/ucsc/bedgraphtobigwig/main"
include { NORMALISATION_FACTORS_REPORT } from "../../modules/local/normalisation_factors_report"
include { NORMALISATION_SCOPE_REFERENCE_REPORT } from "../../modules/local/normalisation_scope_reference_report"

workflow PREPARE_PEAKCALLING {
    take:
    ch_bam         // channel: [ val(meta), [ bam ] ]
    ch_bai         // channel: [ val(meta), [ bai ] ]
    ch_chrom_sizes // channel: [ sizes ]
    ch_dummy_file  // channel: [ dummy ]
    norm_mode      // value:   ["Spikein", "RPKM", "CPM", "BPM", "RPGC", "None" ]
    metadata       // channel  [ csv ]
    normalisation_scope // value: ["all", "group", "group_condition"]
    igg_scale_scope     // value: ["legacy", "group_condition", "sample"]
    ch_flagstat         // channel: [ val(meta), [ flagstat ] ]

    main:
    ch_versions = Channel.empty()
    ch_bedgraph = Channel.empty()
    def norm_scope = normalisation_scope ?: 'all'
    def igg_scope  = igg_scale_scope ?: 'legacy'
    def median = { List values ->
        if (!values || values.size() == 0) {
            return 0
        }
        def sorted = values.collect { it as Double }.sort()
        def size = sorted.size()
        def mid = (int) (size / 2)
        if (size % 2 == 1) {
            return sorted[mid]
        }
        return (sorted[mid - 1] + sorted[mid]) / 2
    }

    if (norm_mode == "Spikein") {
        /*
        * CHANNEL: Load up alignment metadata into channel
        */
        metadata.splitCsv ( header:true, sep:"," )
            .map { row -> [ row[0].id, row[1] ]}
            .set { ch_metadata }
        //ch_metadata | view

        /*
        * CHANNEL: Calculate scale factor for each sample based on spike-in reads.
        */
        ch_bam.map { row -> [ row[0].id, row[0], row[1] ] }
            .join ( ch_metadata )
            .map { row ->
                def meta = row[1]
                def bam = row[2]
                def denominator = row[3].find{ it.key == "bt2_total_aligned" }?.value?.toString()?.toInteger() ?: 0
                def scope_id = norm_scope == 'group' ? meta.group : (norm_scope == 'group_condition' ? meta.group_condition : 'all')
                [ meta, bam, denominator, scope_id ]
            }
            .set { ch_bam_spikein_reads }

        if (norm_scope == 'all') {
            ch_bam_spikein_reads
                .map { meta, bam, reads, scope_id ->
                    def scale = params.normalisation_c / (reads != 0 ? reads : params.normalisation_c)
                    [ meta, bam, scale, reads, scope_id ]
                }
                .set { ch_bam_scale_factor_report }
        } else {
            ch_bam_spikein_reads
                .map { meta, bam, reads, scope_id -> [ scope_id, reads ] }
                .groupTuple(by: [0])
                .map { scope_id, reads_list -> [ scope_id, median(reads_list) ] }
                .set { ch_scope_ref }

            ch_bam_spikein_reads
                .map { meta, bam, reads, scope_id -> [ scope_id, meta, bam, reads ] }
                .join ( ch_scope_ref )
                .map { scope_id, meta, bam, reads, ref ->
                    def scale = (ref == 0 || reads == 0) ? 1 : ref / reads
                    [ meta, bam, scale, reads, scope_id ]
                }
                .set { ch_bam_scale_factor_report }
        }

        ch_bam_scale_factor_report
            .map { meta, bam, scale, reads, scope_id -> [ meta, bam, scale ] }
            .set { ch_bam_scale_factor }
        // EXAMPLE CHANNEL STRUCT: [id, scale_factor]
        //ch_bam_scale_factor | view

        ch_bam_scale_factor_report
            .map { meta, bam, scale, reads, scope_id ->
                [
                    sample_id: meta.id,
                    group: meta.group,
                    condition: meta.condition,
                    replicate: meta.replicate,
                    spikein_reads: reads,
                    scale_factor: scale,
                    scope_id: scope_id
                ]
            }
            .map { record -> [ record.scope_id, record ] }
            .groupTuple(by: [0])
            .map { scope_id, records -> [ scope_id, records ] }
            .set { ch_norm_factors }

        NORMALISATION_FACTORS_REPORT ( ch_norm_factors )
        ch_versions = ch_versions.mix(NORMALISATION_FACTORS_REPORT.out.versions)

        if (params.dump_scale_factors) {
            def ch_scope_reference = Channel.empty()
            if (norm_scope == 'all') {
                ch_scope_reference = Channel.of([scope_id: 'all', reference_reads: params.normalisation_c])
            } else {
                ch_scope_reference = ch_scope_ref
                    .map { scope_id, ref -> [scope_id: scope_id, reference_reads: ref] }
            }

            ch_scope_reference
                .toList()
                .map { records -> records ?: [] }
                .set { ch_scope_reference_records }

            NORMALISATION_SCOPE_REFERENCE_REPORT ( ch_scope_reference_records )
            ch_versions = ch_versions.mix(NORMALISATION_SCOPE_REFERENCE_REPORT.out.versions)
        }
    }
    else if (norm_mode == "None") {
        /*
        * CHANNEL: Assign scale factor of 1
        */
        ch_bam.map { row ->
                [ row[0], row[1], 1 ]
            }
            .set { ch_bam_scale_factor }
        //ch_bam_scale_factor | view
    }

    if (norm_mode == "Spikein" || norm_mode == "None") {
        /*
        * MODULE: Convert bam files to bedgraph
        */
        BEDTOOLS_GENOMECOV (
            ch_bam_scale_factor,
            ch_dummy_file,
            "bedGraph"
        )
        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)
        ch_bedgraph = BEDTOOLS_GENOMECOV.out.genomecov
        //EXAMPLE CHANNEL STRUCT: [META], BEDGRAPH]
        //BEDTOOLS_GENOMECOV.out.genomecov | view

    } else {
        /*
        * CHANNEL: Combine bam and bai files on id
        */
        ch_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }
        .set { ch_bam_bai }
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI]
        //ch_bam_bai | view

        /*
        * CHANNEL: Extract mapped read counts from flagstat
        */
        ch_flagstat
            .map { row ->
                def meta = row[0]
                def flagstat = row[1]
                def mapped_line = flagstat.text.readLines().find { it.contains(' mapped (') }
                def mapped_reads = mapped_line ? mapped_line.tokenize(' ')[0].toInteger() : 0
                [ meta.id, mapped_reads ]
            }
            .set { ch_flagstat_reads }

        /*
        * CHANNEL: Split files based on igg or not
        */
        ch_bam_bai
            .map { row -> [row[0].id, row].flatten() }
            .join ( ch_flagstat_reads )
            .map { row -> [ row[1], row[2], row[3], row[4] ] }
            .branch { it ->
            target:  it[0].is_control == false
            control: it[0].is_control == true
        }
        .set { ch_bam_bai_split }

        /*
        * CHANNEL: Assign scale factor of 1 to target files
        */
        ch_bam_bai_split.target
            .map { row ->
                [ row[0], row[1], row[2], 1 ]
            }
        .set { ch_bam_bai_split_target }
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI, SCALE_FACTOR]
        //ch_bam_bai_split_target | view

        /*
        * CHANNEL: Assign igg scale factor to target files
        */
        if (igg_scope == 'legacy') {
            ch_bam_bai_split.control
                .map { row ->
                    [ row[0], row[1], row[2], params.igg_scale_factor ]
                }
            .set { ch_bam_bai_split_igg }
        }
        else if (igg_scope == 'group_condition') {
            ch_bam_bai_split.control
                .map { row -> [ "${row[0].group}_${row[0].condition}", row ] }
                .groupTuple(by: [0])
                .flatMap { key, entries ->
                    def reads_list = entries.collect { it[3] }
                    def ref = median(reads_list)
                    entries.collect { entry ->
                        def scale = ref / (entry[3] != 0 ? entry[3] : ref)
                        [ entry[0], entry[1], entry[2], scale ]
                    }
                }
                .set { ch_bam_bai_split_igg }
        }
        else {
            def ch_target_reads_map = ch_bam_bai_split.target
                .map { row -> [ "${row[0].control_group}_${row[0].condition}", row[3] ] }
                .toList()
                .map { list ->
                    def map = [:].withDefault { [] }
                    list.each { entry ->
                        def key = entry[0]
                        map[key] = (map[key] ?: []) + [entry[1]]
                    }
                    map
                }

            ch_bam_bai_split.control
                .combine(ch_target_reads_map)
                .map { row, target_map ->
                    def meta = row[0]
                    def reads = row[3]
                    def key = "${meta.group}_${meta.condition}"
                    def target_reads = target_map.get(key, [])
                    if (!target_reads) {
                        target_reads = target_map.findAll { it.key.startsWith("${meta.group}_") }.values().flatten()
                    }
                    def ref = target_reads ? median(target_reads) : reads
                    def scale = ref / (reads != 0 ? reads : ref)
                    [ meta, row[1], row[2], scale ]
                }
                .set { ch_bam_bai_split_igg }
        }
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI, SCALE_FACTOR]
        //ch_bam_bai_split_igg | view

        /*
        * CHANNEL: Mix the split channels back up
        */
        ch_bam_bai_split_target
            .mix(ch_bam_bai_split_igg)
        .set { ch_bam_bai_scale_factor }
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI, SCALE_FACTOR]
        //ch_bam_bai_scale_factor | view

        /*
        * MODULE: Convert bam files to bedgraph and normalise
        */
        DEEPTOOLS_BAMCOVERAGE (
            ch_bam_bai_scale_factor
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions)
        ch_bedgraph = DEEPTOOLS_BAMCOVERAGE.out.bedgraph
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI]
        //ch_bedgraph | view

    }

    /*
    * MODULE: Sort bedgraph
    */
    BEDTOOLS_SORT (
        ch_bedgraph,
        "bedGraph",
        []
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    /*
    * MODULE: Clip off bedgraphs so none overlap beyond chromosome edge
    */
    UCSC_BEDCLIP (
        BEDTOOLS_SORT.out.sorted,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(UCSC_BEDCLIP.out.versions)
    //EXAMPLE CHANNEL STRUCT: [META], BEDGRAPH]
    //UCSC_BEDCLIP.out.bedgraph | view

    /*
    * MODULE: Convert bedgraph to bigwig
    */
    UCSC_BEDGRAPHTOBIGWIG (
        UCSC_BEDCLIP.out.bedgraph,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions)
    //EXAMPLE CHANNEL STRUCT: [[META], BIGWIG]
    //UCSC_BEDGRAPHTOBIGWIG.out.bigwig | view

    emit:
    bedgraph = UCSC_BEDCLIP.out.bedgraph        // channel: [ val(meta), [ bedgraph ] ]
    bigwig   = UCSC_BEDGRAPHTOBIGWIG.out.bigwig // channel: [ val(meta), [ bigwig ] ]
    versions = ch_versions                      // channel: [ versions.yml ]
}
