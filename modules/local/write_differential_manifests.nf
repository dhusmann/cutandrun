process WRITE_DIFFERENTIAL_MANIFESTS {
    tag "differential_manifests"
    label 'process_single'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/00_manifests" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "conda-forge::python=3.11"
    container "quay.io/biocontainers/python:3.8.3"

    input:
    val samples_json
    val peaks_json
    val run_meta_json
    val outdir
    val chrom_sizes

    output:
    path "differential_manifest.samples.tsv", emit: samples
    path "differential_manifest.peaks.tsv"  , emit: peaks
    path "differential_manifest.run_meta.json", emit: run_meta
    path "chrom_sizes.sizes"               , optional: true, emit: chrom_sizes
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """\
    python - <<'PY'
    import json
    import sys
    import shutil
    from pathlib import Path

    samples = json.loads(r'''${samples_json}''')
    peaks = json.loads(r'''${peaks_json}''')
    run_meta = json.loads(r'''${run_meta_json}''')
    outdir = Path(r'''${outdir}''')

    def resolve_published(path_str):
        if not path_str or path_str == 'NA':
            return path_str
        path = Path(path_str)
        if path.exists() and outdir in path.parents:
            return str(path)
        if not outdir.exists():
            return path_str
        matches = list(outdir.rglob(path.name))
        if matches:
            return str(matches[0])
        return path_str

    samples_sorted = sorted(
        samples,
        key=lambda x: (
            str(x.get('group', '')),
            str(x.get('condition', '')),
            int(x.get('replicate', 0)),
            str(x.get('sample_id', '')),
        ),
    )

    peaks_sorted = sorted(
        peaks,
        key=lambda x: (
            str(x.get('group', '')),
            str(x.get('condition', '')),
            int(x.get('replicate', 0)),
            str(x.get('sample_id', '')),
            str(x.get('caller', '')),
        ),
    )

    with open('differential_manifest.samples.tsv', 'w') as handle:
        handle.write("\\t".join([
            'sample_id',
            'group',
            'condition',
            'replicate',
            'control_group',
            'control_condition',
            'final_bam',
            'final_bai',
            'normalisation_mode',
            'spikein_scale_factor',
            'bigwig',
        ]) + "\\n")
        for row in samples_sorted:
            handle.write("\\t".join([
                str(row.get('sample_id', '')),
                str(row.get('group', '')),
                str(row.get('condition', '')),
                str(row.get('replicate', '')),
                str(row.get('control_group', 'NA')),
                str(row.get('control_condition', 'NA')),
                resolve_published(str(row.get('final_bam', ''))),
                resolve_published(str(row.get('final_bai', ''))),
                str(row.get('normalisation_mode', '')),
                str(row.get('spikein_scale_factor', 'NA')),
                resolve_published(str(row.get('bigwig', 'NA'))),
            ]) + "\\n")

    with open('differential_manifest.peaks.tsv', 'w') as handle:
        handle.write("\\t".join([
            'sample_id',
            'caller',
            'peaks_path',
            'peaks_format',
            'caller_role',
        ]) + "\\n")
        for row in peaks_sorted:
            handle.write("\\t".join([
                str(row.get('sample_id', '')),
                str(row.get('caller', '')),
                resolve_published(str(row.get('peaks_path', ''))),
                str(row.get('peaks_format', 'other')),
                str(row.get('caller_role', '')),
            ]) + "\\n")

    with open('differential_manifest.run_meta.json', 'w') as handle:
        json.dump(run_meta, handle, indent=2, sort_keys=True)
        handle.write("\\n")

    chrom = '${chrom_sizes}'
    if chrom and Path(chrom).exists():
        shutil.copy(chrom, 'chrom_sizes.sizes')
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """.stripIndent()
}
