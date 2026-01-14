process MAKE_DIFFBIND_SAMPLESHEET {
    tag "${group}.${caller}"
    label 'process_single'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/01_diffbind/00_samplesheets/${caller}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "conda-forge::python=3.11"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    val samples_json
    val group
    val caller

    output:
    tuple val(group), val(caller), path("*.diffbind.csv"), emit: samplesheet
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${group}.diffbind"
    """\
    python - <<'PY'
    import csv
    import json

    samples = json.loads(r'''${samples_json}''')
    out = f"${prefix}.csv"

    headers = [
        'SampleID','Tissue','Factor','Condition','Replicate','bamReads','Peaks','PeakCaller','SpikeinScaleFactor'
    ]
    with open(out, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        for sample in samples:
            writer.writerow([
                sample['sample_id'],
                'CUTRUN',
                sample['group'],
                sample['condition'],
                sample['replicate'],
                sample['bam'],
                sample['peaks'],
                sample['caller'],
                sample.get('spikein_scale_factor', 'NA'),
            ])
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """.stripIndent()
}
