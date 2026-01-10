process DIFFERENTIAL_SUMMARY_MERGE {
    tag "differential_summary"
    label 'process_single'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/multiqc" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "conda-forge::python=3.11"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path summary_files, stageAs: 'summaries/*'
    val skipped_json
    path summary_header
    path skipped_header

    output:
    path "differential.summary.all_methods.tsv", emit: summary
    path "differential.skipped.tsv"           , emit: skipped
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """\
    mkdir -p summaries
    : > summary_files.list
    find summaries -type f -name "*.tsv" -print | sort > summary_files.list

    python - <<'PY'
    import csv
    import json
    import sys
    from pathlib import Path

    summary_files = []
    with open('summary_files.list') as handle:
        for line in handle:
            line = line.strip()
            if line:
                summary_files.append(Path(line))
    summary_files = sorted(summary_files)

    with open('${summary_header}', 'r') as header_handle:
        header_lines = header_handle.read().rstrip('\\n')

    rows = []
    for path in summary_files:
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter='\\t')
            for row in reader:
                rows.append(row)

    base_fields = [
        'method','group','caller','treated_condition','control_condition','n_tested','n_fdr_pass','n_up','n_down','use_spikein','span_mode_used'
    ]
    if rows:
        fieldnames = list(base_fields)
        for row in rows:
            for key in row.keys():
                if key not in fieldnames:
                    fieldnames.append(key)
    else:
        fieldnames = list(base_fields)

    with open('differential.summary.all_methods.tsv', 'w', newline='') as out_handle:
        if header_lines:
            out_handle.write(header_lines + '\\n')
        writer = csv.DictWriter(out_handle, delimiter='\\t', fieldnames=fieldnames, lineterminator='\\n')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    skipped = json.loads(r'''${skipped_json}''') if '${skipped_json}' else []

    with open('${skipped_header}', 'r') as header_handle:
        skipped_header = header_handle.read().rstrip('\\n')

    with open('differential.skipped.tsv', 'w', newline='') as out_handle:
        if skipped_header:
            out_handle.write(skipped_header + '\\n')
        fieldnames = ['method','group','caller','reason','details']
        writer = csv.DictWriter(out_handle, delimiter='\\t', fieldnames=fieldnames, lineterminator='\\n')
        writer.writeheader()
        for row in skipped:
            writer.writerow({
                'method': row.get('method', 'NA'),
                'group': row.get('group', 'NA'),
                'caller': row.get('caller', 'NA'),
                'reason': row.get('reason', 'NA'),
                'details': row.get('details', '')
            })
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """.stripIndent()
}
