process CHIPBINNER_COUNTS {
    tag "${group}"
    label 'process_high'

    publishDir = [
        path: { "${params.outdir}/03_peak_calling/06_differential/02_chipbinner/${group}/counts" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    conda "bioconda::bedtools=2.31.1 conda-forge::python=3.11 conda-forge::pyyaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    path bins
    val samples_json
    path bams
    path bais
    path control_bams
    path control_bais
    val control_conditions_json
    val use_input
    val group
    val use_spikein
    val pseudocount
    val ms_coeffs

    output:
    path "chipbinner.bin_counts.tsv"        , emit: counts
    path "chipbinner.normalized_matrix.tsv", emit: normalized
    path "chipbinner.normalization.json"   , emit: norm_info
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """\
    printf "%s\n" ${bams} > bam_paths.txt
    printf "%s\n" ${bais} > bai_paths.txt
    printf "%s\n" ${control_bams} > control_bam_paths.txt
    printf "%s\n" ${control_bais} > control_bai_paths.txt

    python - <<'PY'
    import json
    import subprocess
    import sys
    from pathlib import Path

    samples = json.loads(r'''${samples_json}''')
    with open('bam_paths.txt') as handle:
        bam_paths = [line.strip() for line in handle if line.strip()]
    with open('bai_paths.txt') as handle:
        bai_paths = [line.strip() for line in handle if line.strip()]
    sample_ids = [s['sample_id'] for s in samples]
    scale_factors = [s.get('spikein_scale_factor') for s in samples]
    ms_coeffs_path = r'''${ms_coeffs}'''
    control_conditions = json.loads(r'''${control_conditions_json}''')
    use_input = str('${use_input}').lower() in ("true", "1", "yes")
    use_spikein = str('${use_spikein}').lower() in ("true", "1", "yes")

    if len(bam_paths) != len(sample_ids):
        sys.stderr.write(f"Expected {len(sample_ids)} BAMs but found {len(bam_paths)}\\n")
        sys.exit(1)

    def ensure_bai_for_bams(bams, bais, label):
        missing = []
        for bam in bams:
            bam_path = Path(bam)
            candidates = [
                Path(f"{bam}.bai"),
                bam_path.with_suffix(".bai"),
            ]
            if any(candidate.exists() for candidate in candidates):
                continue
            missing.append(str(bam_path))
        if missing:
            sys.stderr.write(f"Missing BAM index for {label}: {', '.join(missing)}\\n")
            sys.stderr.write(f"Available BAI files: {', '.join(bais) if bais else 'none'}\\n")
            sys.exit(1)

    ensure_bai_for_bams(bam_paths, bai_paths, "chipbinner samples")

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

    def parse_ms_coeffs_file(path):
        try:
            import yaml
        except Exception as exc:
            raise ValueError(f"PyYAML is required to parse MS coefficients: {exc}")
        with open(path) as handle:
            data = yaml.safe_load(handle)
        if data is None:
            return {}

        def parse_entries(entries):
            coeffs = {}
            if not isinstance(entries, list):
                raise ValueError("MS coefficients 'samples' must be a list")
            for entry in entries:
                if not isinstance(entry, dict):
                    raise ValueError("MS coefficients entry must be a mapping")
                sample_id = entry.get("sample_id") or entry.get("sample")
                if not sample_id:
                    raise ValueError("MS coefficients entry missing sample_id")
                value = entry.get("ms_coeff") or entry.get("coefficient") or entry.get("scale_factor") or entry.get("scale")
                if value is None or value == "":
                    raise ValueError(f"MS coefficients entry missing coefficient for {sample_id}")
                coeffs[str(sample_id)] = float(value)
            return coeffs

        if isinstance(data, list):
            return parse_entries(data)
        if isinstance(data, dict):
            if "samples" in data:
                samples = data["samples"]
                if isinstance(samples, list):
                    return parse_entries(samples)
                if isinstance(samples, dict):
                    return {str(key): float(value) for key, value in samples.items()}
                raise ValueError("MS coefficients 'samples' must be a list or mapping")
            return {str(key): float(value) for key, value in data.items()}
        raise ValueError("Unsupported MS coefficients YAML structure")

    ms_coeffs = {}
    if ms_coeffs_path:
        try:
            ms_coeffs = parse_ms_coeffs_file(ms_coeffs_path)
        except Exception as exc:
            sys.stderr.write(f"Failed to parse MS coefficients YAML: {exc}\\n")
            sys.exit(1)
        if not ms_coeffs:
            sys.stderr.write("MS coefficients YAML contained no entries\\n")
            sys.exit(1)
        missing = [sid for sid in sample_ids if sid not in ms_coeffs]
        if missing:
            sys.stderr.write(f"MS coefficients missing entries for samples: {', '.join(missing)}\\n")
            sys.exit(1)

    def counts_from_lines(lines):
        counts = []
        for line in lines:
            if not line.strip():
                continue
            parts = line.split("\\t")
            counts.append([float(x) for x in parts[3:]])
        return counts

    counts = counts_from_lines(lines)
    if not counts:
        counts = []

    scaling_source = "library_size"
    scale_values = None
    if ms_coeffs:
        scaling_source = "ms"
        scale_values = [float(ms_coeffs[sid]) for sid in sample_ids]
    elif use_spikein and all(sf not in (None, "NA", "") for sf in scale_factors):
        scaling_source = "spikein"
        scale_values = [float(sf) for sf in scale_factors]
    else:
        if counts:
            lib_sizes = [0.0 for _ in sample_ids]
            for row in counts:
                for idx, val in enumerate(row):
                    lib_sizes[idx] += val
            lib_sorted = sorted(lib_sizes)
            median_lib = lib_sorted[len(lib_sorted) // 2] if lib_sorted else 1.0
            scale_values = [float(median_lib / ls) if ls > 0 else 1.0 for ls in lib_sizes]
        else:
            scale_values = [1.0 for _ in sample_ids]

    scaled_counts = []
    for row in counts:
        scaled_counts.append([val * scale_values[idx] for idx, val in enumerate(row)])

    if use_input and control_conditions:
        with open('control_bam_paths.txt') as handle:
            control_bam_paths = [line.strip() for line in handle if line.strip()]
        with open('control_bai_paths.txt') as handle:
            control_bai_paths = [line.strip() for line in handle if line.strip()]
        ensure_bai_for_bams(control_bam_paths, control_bai_paths, "chipbinner controls")
        if len(control_bam_paths) != len(control_conditions):
            sys.stderr.write("Control BAMs and conditions mismatch\\n")
            sys.exit(1)
        control_cmd = ["bedtools", "multicov", "-bams"] + control_bam_paths + ["-bed", "${bins}"]
        control_res = subprocess.run(control_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if control_res.returncode != 0:
            sys.stderr.write(control_res.stderr)
            sys.exit(control_res.returncode)
        control_lines = control_res.stdout.strip().splitlines()
        control_counts = counts_from_lines(control_lines)
        if control_counts and len(control_counts) != len(scaled_counts):
            sys.stderr.write("Control counts row mismatch with target counts\\n")
            sys.exit(1)

        cond_to_idx = {cond: idx for idx, cond in enumerate(control_conditions)}
        for idx, sample in enumerate(samples):
            cond = sample.get('condition')
            if cond not in cond_to_idx:
                sys.stderr.write(f"Missing control counts for condition {cond}\\n")
                sys.exit(1)
            ctrl_idx = cond_to_idx[cond]
            for row_idx, row in enumerate(scaled_counts):
                row[idx] = row[idx] - (control_counts[row_idx][ctrl_idx] * scale_values[idx])
        for row in scaled_counts:
            for idx, val in enumerate(row):
                row[idx] = val if val > 0 else 0.0

    pseudocount = float('${pseudocount}')
    for row in scaled_counts:
        for idx, val in enumerate(row):
            row[idx] = val + pseudocount

    with open('chipbinner.normalized_matrix.tsv', 'w') as handle:
        handle.write("\\t".join(["chr", "start", "end"] + sample_ids) + "\\n")
        for idx, line in enumerate(lines):
            if not line.strip():
                continue
            parts = line.split("\\t")
            coords = parts[:3]
            row = scaled_counts[idx] if scaled_counts else []
            handle.write("\\t".join(coords + [f"{c:.6f}" for c in row]) + "\\n")

    norm_info = {
        "use_input": bool(use_input and control_conditions),
        "use_ms_scaling": scaling_source == "ms",
        "use_spikein_scaling": scaling_source == "spikein",
        "use_library_size": scaling_source == "library_size",
        "scaling_source": scaling_source,
        "control_conditions": control_conditions if use_input and control_conditions else [],
    }
    with open('chipbinner.normalization.json', 'w') as handle:
        json.dump(norm_info, handle, indent=2)
        handle.write("\\n")
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        python: \$(python --version | awk '{print \$2}')
    END_VERSIONS
    """.stripIndent()
}
