process SPAN_CAPABILITY_PROBE {
    tag "omnipeaks"
    label 'process_single'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/03_span/00_manifests" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "conda-forge::openjdk=21.0.2 conda-forge::python=3.11"
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
    python - <<'PY'
    import json
    import re

    commands = set()
    with open('omnipeaks_help.txt') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            m = re.match(r'^(?:[A-Z]+)?\s*([a-zA-Z][\w-]+)\b', line)
            if m:
                cmd = m.group(1)
                if cmd not in ('usage', 'options'):
                    commands.add(cmd)
    commands = sorted(commands)
    has_compare = any(cmd in ('compare', 'diff', 'differential') for cmd in commands)
    out = {
        "commands": commands,
        "has_compare": has_compare,
    }
    with open('omnipeaks_capabilities.json', 'w') as handle:
        json.dump(out, handle, indent=2)
        handle.write("\\n")
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed -e 's/"//g')
    END_VERSIONS
    """.stripIndent()
}
