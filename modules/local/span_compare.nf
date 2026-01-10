process SPAN_COMPARE {
    tag "${meta.id}"
    label 'process_long'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/03_span/${meta.group}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "conda-forge::openjdk=21.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://eclipse-temurin:21-jre' :
        'eclipse-temurin:21-jre' }"

    input:
    tuple val(meta), path(treated_bams), path(control_bams)
    path chrom_sizes
    path omnipeaks_jar
    val gap
    val bin
    val fdr
    val java_heap
    val compare_cmd
    val multibam_mode

    output:
    tuple val(meta), path("span.differential.tsv")     , emit: tsv
    tuple val(meta), path("span.differential.bed")     , emit: bed
    tuple val(meta), path("span.significant.bed")      , emit: significant
    tuple val(meta), path("span.significant_up.bed")   , emit: up
    tuple val(meta), path("span.significant_down.bed") , emit: down
    tuple val(meta), path("span.summary.tsv")          , emit: summary
    path "00_manifests/span_compare_inputs.tsv"        , emit: manifest
    path "span.readme.txt"                             , optional: true, emit: readme
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.span_compare"
    def fdr_arg = fdr ? "--fdr ${fdr}" : ''
    def bin_arg = bin ? "--bin ${bin}" : ''
    def treated_manifest = meta.treated_bams instanceof List ? meta.treated_bams.join(';') : (meta.treated_bams ?: '')
    def control_manifest = meta.control_bams instanceof List ? meta.control_bams.join(';') : (meta.control_bams ?: '')
    def input_mode = meta.input_mode ?: 'native'
    def compare_mode = multibam_mode ?: 'single'
    def compare_command = compare_cmd ?: 'compare'
    """
    java_cmd=\${JAVA_CMD:-java}
    compare_mode="${compare_mode}"
    compare_command="${compare_command}"

    treated_list=( ${treated_bams} )
    control_list=( ${control_bams} )
    if [ "\${compare_mode}" = "repeat" ]; then
        treated_args=""
        for bam in "\${treated_list[@]}"; do
            treated_args="\${treated_args} -t \${bam}"
        done
        control_args=""
        for bam in "\${control_list[@]}"; do
            control_args="\${control_args} -c \${bam}"
        done
    else
        treated_joined=\$(IFS=,; echo "\${treated_list[*]}")
        control_joined=\$(IFS=,; echo "\${control_list[*]}")
        treated_args="-t \${treated_joined}"
        control_args="-c \${control_joined}"
    fi

    "\$java_cmd" -Xmx${java_heap} -jar ${omnipeaks_jar} \${compare_command} \
        \${treated_args} \
        \${control_args} \
        --cs ${chrom_sizes} \
        --gap ${gap} \
        ${bin_arg} \
        ${fdr_arg} \
        -p ${prefix}

    diff_file=\$(ls ${prefix}*.peak 2>/dev/null | head -n 1 || true)
    if [ -z "\$diff_file" ]; then
        diff_file=\$(ls ${prefix}*.bed 2>/dev/null | head -n 1 || true)
    fi
    if [ -z "\$diff_file" ]; then
        echo "SPAN compare produced no output" >&2
        exit 1
    fi

    cp "\$diff_file" span.differential.tsv

    cat <<'AWK' > span_compare_bed.awk
    function isnum(x) { return (x ~ /^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?\$/) }
    BEGIN { OFS="\\t"; header_done=0; chr_idx=1; start_idx=2; end_idx=3 }
    /^#/ || /^track/ || /^browser/ { next }
    !header_done {
        if (!isnum(\$2)) {
            for (i=1; i<=NF; i++) {
                col=tolower(\$i)
                if (col ~ /^(chr|chrom|seqnames|seqname)\$/) chr_idx=i
                if (col == "start") start_idx=i
                if (col == "end") end_idx=i
            }
            header_done=1
            next
        } else {
            header_done=1
        }
    }
    { print \$chr_idx,\$start_idx,\$end_idx }
    AWK
    awk -f span_compare_bed.awk "\$diff_file" > span.differential.bed

    : > span.significant_up.bed
    : > span.significant_down.bed
    cat <<'AWK' > span_compare_filter.awk
    function isnum(x) { return (x ~ /^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?\$/) }
    function calc_fdr(val) {
        if (!isnum(val)) return 1
        if (val <= 1) return val + 0
        return 10^(-val)
    }
    function calc_log2(val) {
        if (!isnum(val)) return ""
        if (val < 0) return val + 0
        if (val == 0) return 0
        return log(val) / log(2)
    }
    BEGIN { OFS="\\t"; header_done=0; chr_idx=1; start_idx=2; end_idx=3; log2_idx=0; fdr_idx=0; qval_idx=0; fold_idx=0; proxy_log2=0 }
    /^#/ || /^track/ || /^browser/ { next }
    !header_done {
        if (!isnum(\$2)) {
            for (i=1; i<=NF; i++) {
                col=tolower(\$i)
                if (col ~ /^(chr|chrom|seqnames|seqname)\$/) chr_idx=i
                if (col == "start") start_idx=i
                if (col == "end") end_idx=i
                if (col ~ /log2fc|log2foldchange|log2_fold_change|log2fold|log_fc|log2ratio|log2_ratio/) log2_idx=i
                if (col ~ /^(fdr|false.discovery|padj|adj_p|adjp)\$/) fdr_idx=i
                if (col ~ /qvalue|qval|q-value/) qval_idx=i
                if ((col ~ /fold|fc/) && log2_idx==0) fold_idx=i
            }
            header_done=1
            next
        } else {
            chr_idx=1; start_idx=2; end_idx=3; log2_idx=7; qval_idx=9
            proxy_log2=1
            header_done=1
        }
    }
    {
        log2_val=""
        if (log2_idx>0 && isnum(\$log2_idx)) log2_val=\$log2_idx+0
        if (log2_val=="" && fold_idx>0 && isnum(\$fold_idx)) { log2_val=calc_log2(\$fold_idx); proxy_log2=1 }
        if (log2_val=="" && log2_idx==7 && isnum(\$7)) { log2_val=calc_log2(\$7); proxy_log2=1 }
        fdr_val=""
        if (fdr_idx>0) fdr_val=\$fdr_idx
        else if (qval_idx>0) fdr_val=\$qval_idx
        else if (NF>=9) fdr_val=\$9
        fdr_val=calc_fdr(fdr_val)
        if (fdr_val <= fdr_thresh) {
            if (log2_val > 0) print \$chr_idx,\$start_idx,\$end_idx >> up_file
            else if (log2_val < 0) print \$chr_idx,\$start_idx,\$end_idx >> down_file
        }
    }
    END {
        if (proxy_log2) print "proxy_log2" > proxy_file
    }
    AWK

    awk -v fdr_thresh="${fdr}" -v up_file="span.significant_up.bed" -v down_file="span.significant_down.bed" -v proxy_file="span.proxy_log2.txt" -f span_compare_filter.awk "\$diff_file"

    cat span.significant_up.bed span.significant_down.bed | awk 'NF' | sort -k1,1 -k2,2n | uniq > span.significant.bed

    total=\$(awk 'function isnum(x) { return (x ~ /^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?\$/) } /^#/ || /^track/ || /^browser/ { next } !header_seen { if (!isnum(\$2)) { header_seen=1; next } header_seen=1 } { n++ } END { print n+0 }' "\$diff_file")
    n_up=\$(wc -l < span.significant_up.bed | awk '{print \$1}')
    n_down=\$(wc -l < span.significant_down.bed | awk '{print \$1}')
    n_fdr_pass=\$((n_up + n_down))
    cat <<-END_SUMMARY > span.summary.tsv
    method\tgroup\tcaller\ttreated_condition\tcontrol_condition\tn_tested\tn_fdr_pass\tn_up\tn_down\tuse_spikein\tspan_mode_used
    span\t${meta.group}\tNA\t${meta.treated ?: 'NA'}\t${meta.control ?: 'NA'}\t\${total}\t\${n_fdr_pass}\t\${n_up}\t\${n_down}\tNA\tnative
    END_SUMMARY

    mkdir -p 00_manifests
    if [ -z "${treated_manifest}" ]; then
        treated_manifest=\$(IFS=';'; echo "\${treated_list[*]}")
    else
        treated_manifest="${treated_manifest}"
    fi
    if [ -z "${control_manifest}" ]; then
        control_manifest=\$(IFS=';'; echo "\${control_list[*]}")
    else
        control_manifest="${control_manifest}"
    fi
    {
        printf "group\\tcondition\\tinput_bams\\tinput_mode\\n"
        printf "%s\\t%s\\t%s\\t%s\\n" "${meta.group}" "${meta.treated ?: 'NA'}" "\${treated_manifest}" "${input_mode}"
        printf "%s\\t%s\\t%s\\t%s\\n" "${meta.group}" "${meta.control ?: 'NA'}" "\${control_manifest}" "${input_mode}"
    } > 00_manifests/span_compare_inputs.tsv

    if [ -f span.proxy_log2.txt ]; then
        cat <<-END_README > span.readme.txt
    SPAN native differential parsing notes:
    - No explicit log2FC column detected; log2FC was derived from the signal/fold-change column.
    - For positive values, log2FC = log2(value); negative values are treated as already log2-scaled.
    - FDR was computed from the Q-value column (column 9) assuming -log10 scale when values > 1.
    - Direction follows treated (${meta.treated ?: 'NA'}) vs control (${meta.control ?: 'NA'}).
    END_README
    fi
    rm -f span_compare_filter.awk span_compare_bed.awk span.proxy_log2.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(\${java_cmd} -version 2>&1 | head -n 1 | sed -e 's/\"//g')
    END_VERSIONS
    """
}
