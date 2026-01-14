process SPAN_CAPABILITY_PROBE {
    tag "omnipeaks"
    label 'process_single'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/03_span/00_manifests" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "conda-forge::openjdk=21.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://eclipse-temurin:21-jre' :
        'eclipse-temurin:21-jre' }"

    input:
    path omnipeaks_jar

    output:
    path "omnipeaks_capabilities.json", emit: capabilities
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """\
    java -jar ${omnipeaks_jar} --help > omnipeaks_help.txt
    commands=\$(awk '
        {
            line=\$0
            sub(/^[ \\t]+/, "", line)
            if (line == "") next
            if (match(line, /^([A-Za-z][A-Za-z0-9_-]*)[ \\t]+/, m)) {
                cmd=m[1]
                lc=tolower(cmd)
                if (lc != "usage" && lc != "options" && lc != "commands") print cmd
            }
        }' omnipeaks_help.txt | sort -u)

    has_compare=false
    compare_cmd=""
    for cmd in \${commands}; do
        lc=\$(printf '%s' "\${cmd}" | tr '[:upper:]' '[:lower:]')
        if [ "\${lc}" = "compare" ] || [ "\${lc}" = "diff" ] || [ "\${lc}" = "differential" ]; then
            has_compare=true
        fi
        if [ -z "\${compare_cmd}" ] && [ "\${lc}" = "compare" ]; then
            compare_cmd="\${cmd}"
        fi
    done
    if [ -z "\${compare_cmd}" ]; then
        for cmd in \${commands}; do
            lc=\$(printf '%s' "\${cmd}" | tr '[:upper:]' '[:lower:]')
            if [ "\${lc}" = "diff" ] || [ "\${lc}" = "differential" ]; then
                compare_cmd="\${cmd}"
                break
            fi
        done
    fi

    compare_multibam=false
    compare_multibam_mode="single"
    if [ -n "\${compare_cmd}" ]; then
        java -jar ${omnipeaks_jar} "\${compare_cmd}" --help > omnipeaks_compare_help.txt 2>/dev/null || true
        if [ -s omnipeaks_compare_help.txt ]; then
            if grep -Eiq 'comma[- ]?separated|comma separated|list|multiple|replicate|replicates|bams?' omnipeaks_compare_help.txt; then
                compare_multibam=true
            fi
            if grep -Eiq 'comma[- ]?separated|bam[^ ]*,|bam\[,?bam' omnipeaks_compare_help.txt; then
                compare_multibam=true
                compare_multibam_mode="comma"
            fi
            if grep -Eiq 'repeat|multiple[[:space:]]+times|use[[:space:]]+multiple' omnipeaks_compare_help.txt; then
                if [ "\${compare_multibam_mode}" = "single" ]; then
                    compare_multibam_mode="repeat"
                fi
                compare_multibam=true
            fi
            if [ "\${compare_multibam}" = true ] && [ "\${compare_multibam_mode}" = "single" ]; then
                compare_multibam_mode="comma"
            fi
        fi
    fi

    if [ -n "\${commands}" ]; then
        commands_json=\$(printf '%s\n' "\${commands}" | awk 'BEGIN{first=1} {gsub(/\"/,"\\\\\""); if (!first) printf ", "; printf "\"%s\"", \$0; first=0} END{print ""}')
        commands_json="[ \${commands_json} ]"
    else
        commands_json="[]"
    fi

    cat <<-END_JSON > omnipeaks_capabilities.json
    {
      "commands": \${commands_json},
      "has_compare": \${has_compare},
      "compare_command": "\${compare_cmd}",
      "compare_multi_bam": \${compare_multibam},
      "compare_multi_bam_mode": "\${compare_multibam_mode}"
    }
    END_JSON

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed -e 's/"//g')
    END_VERSIONS
    """.stripIndent()
}
