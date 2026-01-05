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
    for cmd in \${commands}; do
        lc=\$(printf '%s' "\${cmd}" | tr '[:upper:]' '[:lower:]')
        if [ "\${lc}" = "compare" ] || [ "\${lc}" = "diff" ] || [ "\${lc}" = "differential" ]; then
            has_compare=true
            break
        fi
    done

    if [ -n "\${commands}" ]; then
        commands_json=\$(printf '%s\n' "\${commands}" | awk 'BEGIN{first=1} {gsub(/\"/,"\\\\\""); if (!first) printf ", "; printf "\"%s\"", \$0; first=0} END{print ""}')
        commands_json="[ \${commands_json} ]"
    else
        commands_json="[]"
    fi

    cat <<-END_JSON > omnipeaks_capabilities.json
    {
      "commands": \${commands_json},
      "has_compare": \${has_compare}
    }
    END_JSON

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed -e 's/"//g')
    END_VERSIONS
    """.stripIndent()
}
