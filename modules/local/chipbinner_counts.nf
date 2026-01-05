process CHIPBINNER_COUNTS {
    tag "${group}"
    label 'process_high'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/02_chipbinner/${group}/counts" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "bioconda::bedtools=2.31.1 conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    path bins
    val samples_json
    path bams
    path bais
    val group
    val use_spikein
    val pseudocount

    output:
    path "chipbinner.bin_counts.tsv"        , emit: counts
    path "chipbinner.normalized_matrix.tsv", emit: normalized
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """\
    printf "%s\n" ${bams} > bam_paths.txt

    python - <<'PY'
    import json
    import subprocess
    import sys
    from pathlib import Path

    samples = json.loads(r'''${samples_json}''')
    with open('bam_paths.txt') as handle:
        bam_paths = [line.strip() for line in handle if line.strip()]
    sample_ids = [s['sample_id'] for s in samples]
    scale_factors = [s.get('spikein_scale_factor') for s in samples]
    if len(bam_paths) != len(sample_ids):
        sys.stderr.write(f\"Expected {len(sample_ids)} BAMs but found {len(bam_paths)}\\n\")
        sys.exit(1)

    cmd = ["bedtools", "multicov", "-bams"] + bam_paths + ["-bed", "${bins}"]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        sys.exit(result.returncode)

    lines = result.stdout.strip().splitlines()
    with open('chipbinner.bin_counts.tsv', 'w') as handle:
        handle.write("\\t".join(["chr", "start", "end"] + sample_ids) + "\\n")
        for line in lines:
            if not line.strip():
                continue
            handle.write(line + "\\n")

    # normalized matrix
    use_spikein = str('${use_spikein}').lower() in ("true", "1", "yes")
    with open('chipbinner.normalized_matrix.tsv', 'w') as handle:
        handle.write("\\t".join(["chr", "start", "end"] + sample_ids) + "\\n")
        for line in lines:
            if not line.strip():
                continue
            parts = line.split("\\t")
            coords = parts[:3]
            counts = [float(x) for x in parts[3:]]
            if use_spikein and all(sf not in (None, "NA", "") for sf in scale_factors):
                factors = [float(sf) for sf in scale_factors]
                counts = [c * f for c, f in zip(counts, factors)]
            counts = [c + float('${pseudocount}') for c in counts]
            handle.write("\\t".join(coords + [f"{c:.6f}" for c in counts]) + "\\n")
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """.stripIndent()
}
