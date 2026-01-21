process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    conda "bioconda::multiqc=1.19"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.19--pyhdfd78af_0' :
        'biocontainers/multiqc:1.19--pyhdfd78af_0' }"

    input:
    path versions

    output:
    path "software_versions.yml"           , emit: yml
    path "software_versions_mqc.yml"       , emit: mqc_yml
    path "software_versions_unique_mqc.yml", emit: mqc_unique_yml
    path "local_versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python

    import yaml
    import platform
    import json
    from textwrap import dedent

    def _make_versions_html(versions):
        html = [
            dedent(
                '''\\
                <style>
                #nf-core-versions tbody:nth-child(even) {
                    background-color: #f2f2f2;
                }
                </style>
                <table class="table" style="width:100%" id="nf-core-versions">
                    <thead>
                        <tr>
                            <th> Process Name </th>
                            <th> Software </th>
                            <th> Version  </th>
                        </tr>
                    </thead>
                '''
            )
        ]
        for process, tmp_versions in sorted(versions.items()):
            html.append("<tbody>")
            for i, (tool, version) in enumerate(sorted(tmp_versions.items())):
                html.append(
                    dedent(
                        f'''\\
                        <tr>
                            <td><samp>{process if (i == 0) else ''}</samp></td>
                            <td><samp>{tool}</samp></td>
                            <td><samp>{version}</samp></td>
                        </tr>
                        '''
                    )
                )
            html.append("</tbody>")
        html.append("</table>")
        return "\\n".join(html)

    def _make_versions_unique_html(versions):
        unique_versions = []

        for process, tmp_versions in sorted(versions.items()):
            for i, (tool, version) in enumerate(sorted(tmp_versions.items())):
                tool_version = tool + "=" + version
                if tool_version not in unique_versions:
                    unique_versions.append(tool_version)

        unique_versions.sort()

        html = [
            dedent(
                '''\\
                <style>
                #nf-core-versions-unique tbody:nth-child(even) {
                    background-color: #f2f2f2;
                }
                </style>
                <table class="table" style="width:100%" id="nf-core-versions-unique">
                    <thead>
                        <tr>
                            <th> Software </th>
                            <th> Version  </th>
                        </tr>
                    </thead>
                '''
            )
        ]

        for tool_version in unique_versions:
            tool_version_split = tool_version.split('=')
            html.append("<tbody>")
            html.append(
                dedent(
                    f'''\\
                    <tr>
                        <td><samp>{tool_version_split[0]}</samp></td>
                        <td><samp>{tool_version_split[1]}</samp></td>
                    </tr>
                    '''
                )
            )
            html.append("</tbody>")
        html.append("</table>")
        return "\\n".join(html)

    module_versions = {}
    module_versions["${task.process}"] = {
        'python': platform.python_version(),
        'yaml': yaml.__version__
    }

    with open("$versions") as f:
        raw_lines = f.read().splitlines()

    # Sanitize potentially noisy version strings (eg. warnings) into valid YAML
    clean_lines = []
    for line in raw_lines:
        if not line.strip():
            continue
        if line.startswith('"'):
            clean_lines.append(line)
            continue
        stripped = line.lstrip(' ')
        indent = line[:len(line) - len(stripped)]
        if ':' not in stripped:
            # Drop stray lines (eg. warning continuations)
            continue
        key_part, value = stripped.split(':', 1)
        key = f"{indent}{key_part}:"
        value = value.strip()
        if value == "":
            clean_lines.append(line)
            continue
        if not (value.startswith('"') or value.startswith("'")):
            clean_lines.append(f"{key} {json.dumps(value)}")
        else:
            clean_lines.append(line)

    workflow_versions = yaml.load("\\n".join(clean_lines) + "\\n", Loader=yaml.BaseLoader) | module_versions

    workflow_versions["Workflow"] = {
        "Nextflow": "$workflow.nextflow.version",
        "$workflow.manifest.name": "$workflow.manifest.version"
    }

    versions_mqc = {
        'parent_id': 'software_versions',
        'parent_name': 'Software Versions',
        'parent_description': 'Details software versions used in the pipeline run',
        'id': 'software-versions-by-process',
        'section_name': '${workflow.manifest.name} software versions by process',
        'section_href': 'https://github.com/${workflow.manifest.name}',
        'plot_type': 'html',
        'description': 'are collected at run time from the software output.',
        'data': _make_versions_html(workflow_versions)
    }

    versions_mqc_unique = {
        'parent_id': 'software_versions',
        'parent_name': 'Software Versions',
        'parent_description': 'Details software versions used in the pipeline run',
        'id': 'software-versions-unique',
        'section_name': '${workflow.manifest.name} Software Versions',
        'section_href': 'https://github.com/${workflow.manifest.name}',
        'plot_type': 'html',
        'description': 'are collected at run time from the software output.',
        'data': _make_versions_unique_html(workflow_versions)
    }

    with open("software_versions.yml", 'w') as f:
        yaml.dump(workflow_versions, f, default_flow_style=False)

    with open("software_versions_mqc.yml", 'w') as f:
        yaml.dump(versions_mqc, f, default_flow_style=False)

    with open("software_versions_unique_mqc.yml", 'w') as f:
        yaml.dump(versions_mqc_unique, f, default_flow_style=False)

    with open('local_versions.yml', 'w') as f:
        yaml.dump(module_versions, f, default_flow_style=False)
    """
}
