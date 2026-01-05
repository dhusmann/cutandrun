#!/usr/bin/env python3
import argparse
import csv
import math
import os
import re
import subprocess
import sys
from typing import Dict, List, Tuple


def read_tsv(path: str) -> List[Dict[str, str]]:
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for row in reader]


def write_tsv(path: str, header: List[str], rows: List[Dict[str, str]]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_list(value: str) -> List[str]:
    if not value:
        return []
    return [x.strip().lower() for x in value.split(",") if x.strip()]


def is_span_caller(name: str) -> bool:
    if not name:
        return False
    lowered = str(name).lower()
    return lowered.startswith("span") or lowered.startswith("omnipeak")


def safe_float(value: str):
    if value is None:
        return None
    try:
        return float(value)
    except Exception:
        return None


def run_cmd(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, check=check)


def load_help_text(jar_path: str, java_heap: str, args: List[str]) -> Tuple[bool, str]:
    cmd = ["java", f"-Xmx{java_heap}", "-jar", jar_path] + args
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except FileNotFoundError:
        return False, "java_missing"
    output = (result.stdout or "") + (result.stderr or "")
    return True, output


def detect_signature(jar_path: str, java_heap: str) -> Tuple[bool, str, str]:
    ok, help_text = load_help_text(jar_path, java_heap, ["--help"])
    if not ok:
        return False, "java_missing", ""
    if not help_text:
        return False, "help_failed", ""

    lower_help = help_text.lower()
    supports_compare = False
    if re.search(r"\bcompare\b", lower_help):
        supports_compare = True
    elif re.search(r"\bdifferential\b", lower_help):
        supports_compare = True

    compare_help = ""
    if supports_compare:
        ok_cmp, cmp_text = load_help_text(jar_path, java_heap, ["compare", "--help"])
        if ok_cmp and cmp_text:
            compare_help = cmp_text
            lower_cmp = cmp_text.lower()
            if re.search(r"unknown command|no such command", lower_cmp):
                supports_compare = False
        else:
            compare_help = help_text

    if not supports_compare:
        return False, "no_compare", help_text

    signature = "compare:pooled"
    probe_text = compare_help or help_text
    for line in probe_text.splitlines():
        lline = line.lower()
        if ("treat" in lline or "control" in lline or "case" in lline) and (
            "replicate" in lline
            or "replicates" in lline
            or "comma" in lline
            or "list" in lline
            or "bams" in lline
            or "bam" in lline
            or "files" in lline
        ):
            signature = "compare:list"
            break
    if "replicate" in probe_text.lower() or "comma" in probe_text.lower():
        signature = "compare:list"

    return True, signature, probe_text


def parse_flags(help_text: str) -> Dict[str, str]:
    if not help_text:
        return {}
    long_flags = set(re.findall(r"--[A-Za-z0-9][A-Za-z0-9_.-]*", help_text))
    short_flags = set(re.findall(r"(?<!\w)-[A-Za-z](?!\w)", help_text))

    def pick(candidates: List[str]) -> str:
        for cand in candidates:
            if cand in long_flags or cand in short_flags:
                return cand
        return ""

    flags = {
        "treatment": pick(["--treatment", "--treated", "--treat", "--case", "-t"]),
        "control": pick(["--control", "--ctrl", "--background", "--input", "-c"]),
        "chrom": pick(["--chrom.sizes", "--chrom_sizes", "--chrom", "--cs"]),
        "bin": pick(["--bin", "-b"]),
        "gap": pick(["--gap", "-g"]),
        "fdr": pick(["--fdr", "--q", "--qvalue"]),
        "output": pick(["--peaks", "--output", "--out", "-p", "--prefix"]),
    }
    flags["output_is_prefix"] = ""
    if "--prefix" in long_flags:
        flags["output_is_prefix"] = "true"
    if flags["output"] in {"--prefix"}:
        flags["output_is_prefix"] = "true"
    if flags["output"] == "-p":
        for line in help_text.splitlines():
            if "-p" in line and "prefix" in line.lower():
                flags["output_is_prefix"] = "true"
                break
    if flags["output"] == "--peaks":
        flags["output_is_prefix"] = "false"
    return flags


def is_int(value: str) -> bool:
    try:
        int(value)
        return True
    except Exception:
        return False


def read_intervals(path: str) -> List[Tuple[str, int, int]]:
    intervals = []
    with open(path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            if not is_int(parts[1]) or not is_int(parts[2]):
                continue
            start = int(parts[1])
            end = int(parts[2])
            if end <= start:
                continue
            intervals.append((parts[0], start, end))
    return intervals


def merge_intervals(intervals: List[Tuple[str, int, int]]) -> List[Tuple[str, int, int]]:
    if not intervals:
        return []
    intervals.sort(key=lambda x: (x[0], x[1], x[2]))
    merged = [list(intervals[0])]
    for chrom, start, end in intervals[1:]:
        last = merged[-1]
        if chrom != last[0] or start > last[2]:
            merged.append([chrom, start, end])
        else:
            last[2] = max(last[2], end)
    return [tuple(x) for x in merged]


def write_bed(path: str, intervals: List[Tuple[str, int, int]]) -> None:
    with open(path, "w") as handle:
        for chrom, start, end in intervals:
            handle.write(f"{chrom}\t{start}\t{end}\n")


def parse_native_output(path: str) -> List[Dict[str, str]]:
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path, "r") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip() and not line.startswith("#")]
    if not lines:
        return rows

    header = None
    first = lines[0].split("\t")
    lower = [col.lower() for col in first]
    if any(col in ("chr", "chrom", "chromosome") for col in lower) and any(
        col.startswith("start") or col == "chromstart" for col in lower
    ):
        header = lower
        data_lines = lines[1:]
    else:
        data_lines = lines

    def idx_for(names: List[str]):
        if not header:
            return None
        for name in names:
            if name in header:
                return header.index(name)
        return None

    idx_chr = idx_for(["chr", "chrom", "chromosome", "seqname"]) if header else 0
    idx_start = idx_for(["start", "chromstart", "begin"]) if header else 1
    idx_end = idx_for(["end", "chromend", "stop"]) if header else 2
    idx_log = idx_for(["log2fc", "logfc", "log2foldchange", "log2_fold_change", "log2foldchange"])
    idx_pval = idx_for(["pval", "pvalue", "p-value", "p.value"]) if header else None
    idx_fdr = idx_for(["fdr", "qvalue", "q-value", "padj", "adj.p.val"]) if header else None

    for line in data_lines:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        chr_val = parts[idx_chr] if idx_chr is not None and idx_chr < len(parts) else parts[0]
        start_val = parts[idx_start] if idx_start is not None and idx_start < len(parts) else parts[1]
        end_val = parts[idx_end] if idx_end is not None and idx_end < len(parts) else parts[2]

        log_val = "NA"
        pval_val = "NA"
        fdr_val = "NA"

        if header and idx_log is not None and idx_log < len(parts):
            log_val = parts[idx_log]
        if header and idx_pval is not None and idx_pval < len(parts):
            pval_val = parts[idx_pval]
        if header and idx_fdr is not None and idx_fdr < len(parts):
            fdr_val = parts[idx_fdr]

        if not header:
            numeric = []
            for idx in range(3, len(parts)):
                val = safe_float(parts[idx])
                if val is not None:
                    numeric.append((idx, val))
            for idx, val in numeric:
                if log_val == "NA" and not (0.0 <= val <= 1.0):
                    log_val = parts[idx]
                    continue
                if pval_val == "NA" and 0.0 <= val <= 1.0:
                    pval_val = parts[idx]
                    continue
                if fdr_val == "NA" and 0.0 <= val <= 1.0:
                    fdr_val = parts[idx]
                    continue
            if log_val == "NA" and numeric:
                log_val = parts[numeric[0][0]]

        rows.append({
            "chr": chr_val,
            "start": start_val,
            "end": end_val,
            "log2FC": log_val,
            "pval": pval_val,
            "FDR": fdr_val,
        })
    return rows


def write_summary(path: str, group: str, treated: str, control: str, n_tested: int, n_fdr: int, n_up: int, n_down: int,
                  mode: str, status: str, reason: str) -> None:
    header = [
        "group",
        "treated",
        "control",
        "n_tested",
        "n_fdr_pass",
        "n_up",
        "n_down",
        "mode",
        "status",
        "reason",
    ]
    rows = [{
        "group": group,
        "treated": treated,
        "control": control,
        "n_tested": str(n_tested),
        "n_fdr_pass": str(n_fdr),
        "n_up": str(n_up),
        "n_down": str(n_down),
        "mode": mode,
        "status": status,
        "reason": reason,
    }]
    write_tsv(path, header, rows)


def write_mode(path: str, mode: str, signature: str, notes: List[str] = None) -> None:
    with open(path, "w") as handle:
        handle.write(f"{mode}:{signature}\n")
        if notes:
            for note in notes:
                handle.write(f"{note}\n")


def filter_bed_by_stats(rows: List[Dict[str, str]], fdr_cutoff: float, direction: str) -> List[Tuple[str, int, int]]:
    out = []
    for row in rows:
        fdr = safe_float(row.get("FDR"))
        log2fc = safe_float(row.get("log2FC"))
        if fdr is None or log2fc is None:
            continue
        if fdr > fdr_cutoff:
            continue
        if direction == "up" and log2fc <= 0:
            continue
        if direction == "down" and log2fc >= 0:
            continue
        if not is_int(row.get("start", "")) or not is_int(row.get("end", "")):
            continue
        out.append((row["chr"], int(row["start"]), int(row["end"])))
    return out


def write_failure_outputs(
    outdir: str,
    group: str,
    treated: str,
    control: str,
    mode: str,
    signature: str,
    reason: str,
    allow_partial: bool,
) -> None:
    diff_path = os.path.join(outdir, "span.differential.tsv")
    bed_path = os.path.join(outdir, "span.differential.peaks.bed")
    up_path = os.path.join(outdir, "span.up.bed")
    down_path = os.path.join(outdir, "span.down.bed")
    summary_path = os.path.join(outdir, "span.summary.tsv")
    mode_path = os.path.join(outdir, "span.mode.txt")

    with open(diff_path, "w") as handle:
        handle.write("chr\tstart\tend\tlog2FC\tpval\tFDR\n")
    open(bed_path, "w").close()
    open(up_path, "w").close()
    open(down_path, "w").close()
    status = "SKIP" if allow_partial else "FAIL"
    write_summary(summary_path, group, treated, control, 0, 0, 0, 0, mode, status, reason)
    write_mode(mode_path, mode, signature)


def compute_orientation_log2fc(
    bed_path: str,
    treated_bams: List[str],
    control_bams: List[str],
    pseudocount: float,
    treated_size_factors: List[float] = None,
    control_size_factors: List[float] = None,
) -> Dict[Tuple[str, str, str], List[float]]:
    if not treated_bams or not control_bams:
        raise RuntimeError("orientation_missing_bams")
    cmd = ["bedtools", "multicov", "-bed", bed_path, "-bams"] + treated_bams + control_bams
    result = run_cmd(cmd, check=False)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or "bedtools_multicov_failed")

    def apply_size_factors(counts: List[float], factors: List[float]) -> List[float]:
        if not factors or len(factors) != len(counts):
            return counts
        adjusted = []
        for count, factor in zip(counts, factors):
            if factor is None or factor == 0:
                factor = 1.0
            adjusted.append(count / factor)
        return adjusted

    orientation = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        key = (parts[0], parts[1], parts[2])
        counts = [safe_float(val) or 0.0 for val in parts[3:]]
        t_counts = counts[: len(treated_bams)]
        c_counts = counts[len(treated_bams) :]
        t_counts = apply_size_factors(t_counts, treated_size_factors)
        c_counts = apply_size_factors(c_counts, control_size_factors)
        if not t_counts or not c_counts:
            continue
        t_mean = sum(t_counts) / len(t_counts)
        c_mean = sum(c_counts) / len(c_counts)
        log2fc = math.log2((t_mean + pseudocount) / (c_mean + pseudocount))
        orientation.setdefault(key, []).append(log2fc)
    return orientation


def compute_spikein_factors(samples: List[Dict[str, str]], use_spikein: bool):
    scale_values = []
    size_factors = []
    numeric_present = False
    for row in samples:
        raw = row.get("spikein_scale_factor", "NA")
        scale_values.append(raw)
        numeric = None
        try:
            numeric = float(raw)
            numeric_present = True
        except Exception:
            numeric = None

        if use_spikein:
            if numeric is None or numeric == 0:
                size = 1.0
            else:
                size = 1.0 / numeric
        else:
            size = None
        size_factors.append(size)

    use_spikein_factors = use_spikein and numeric_present
    return scale_values, size_factors, use_spikein_factors


def write_normalization_factors(
    path: str,
    samples: List[Dict[str, str]],
    scale_values: List[str],
    size_factors: List[float],
    use_spikein_factors: bool,
    backend: str,
) -> None:
    with open(path, "w") as handle:
        handle.write("sample_id\tgroup\tcondition\tspikein_scale_factor\tspikein_size_factor\tused_by_backend\n")
        for row, scale_val, size_val in zip(samples, scale_values, size_factors):
            used_by_backend = "spikein" if use_spikein_factors else backend
            size_str = "NA"
            if use_spikein_factors:
                size = size_val if size_val is not None else 1.0
                size_str = f"{size:.6f}"
            handle.write(
                f"{row.get('sample_id','')}\t{row.get('group','')}\t{row.get('condition','')}\t{scale_val}\t{size_str}\t{used_by_backend}\n"
            )


def main():
    parser = argparse.ArgumentParser(description="SPAN/OmniPeak differential runner")
    parser.add_argument("--jar", required=True)
    parser.add_argument("--mode", required=True)
    parser.add_argument("--contrast", required=True)
    parser.add_argument("--group", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--peaks", required=True)
    parser.add_argument("--caller-priority", default="")
    parser.add_argument("--chrom-sizes")
    parser.add_argument("--bin", type=int, default=200)
    parser.add_argument("--gap", type=int, default=5)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--fallback-backend", default="DESeq2")
    parser.add_argument("--use-spikein", default="false")
    parser.add_argument("--allow-partial", default="false")
    parser.add_argument("--java-heap", default="8G")
    parser.add_argument("--cpus", type=int, default=1)
    parser.add_argument("--pooling-dir", default="")
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    contrast = [x.strip() for x in args.contrast.split(",") if x.strip()]
    treated = contrast[0] if len(contrast) > 0 else "treated"
    control = contrast[1] if len(contrast) > 1 else "control"

    use_spikein = str(args.use_spikein).lower() in {"true", "1", "yes"}
    allow_partial = str(args.allow_partial).lower() in {"true", "1", "yes"}

    native_supported, signature, help_text = detect_signature(args.jar, args.java_heap)

    chosen_mode = "fallback"
    if args.mode == "native":
        if not native_supported:
            reason = "native_not_supported"
            write_failure_outputs(args.outdir, args.group, treated, control, "native", signature, reason, allow_partial)
            if allow_partial:
                sys.exit(0)
            print("ERROR: SPAN jar does not support native compare.", file=sys.stderr)
            sys.exit(1)
        chosen_mode = "native"
    elif args.mode == "fallback":
        chosen_mode = "fallback"
    else:
        chosen_mode = "native" if native_supported else "fallback"

    samples = read_tsv(args.samples)
    peaks = read_tsv(args.peaks)

    samples_group = [row for row in samples if row.get("group") == args.group and row.get("condition") in {treated, control}]
    if not samples_group:
        reason = "no_samples"
        write_failure_outputs(args.outdir, args.group, treated, control, chosen_mode, signature, reason, allow_partial)
        if allow_partial:
            sys.exit(0)
        print("ERROR: No samples found for group.", file=sys.stderr)
        sys.exit(1)

    peaks_group = [row for row in peaks if row.get("group") == args.group and row.get("condition") in {treated, control}]
    span_peaks = [row for row in peaks_group if is_span_caller(row.get("caller"))]
    if not span_peaks:
        reason = "no_span_peaks"
        write_failure_outputs(args.outdir, args.group, treated, control, chosen_mode, signature, reason, allow_partial)
        if allow_partial:
            sys.exit(0)
        print("ERROR: No SPAN/OmniPeak peaks found for group.", file=sys.stderr)
        sys.exit(1)

    priority = [x for x in parse_list(args.caller_priority) if is_span_caller(x)]
    available_callers = sorted({(row.get("caller") or "").lower() for row in span_peaks})
    chosen_caller = None
    for cand in priority:
        if cand in available_callers:
            chosen_caller = cand
            break
    if not chosen_caller:
        chosen_caller = available_callers[0]
    if len(available_callers) > 1:
        print(
            f"WARNING: Multiple SPAN callers for group {args.group}; using '{chosen_caller}'.",
            file=sys.stderr,
        )

    peaks_by_sample = {}
    for row in sorted(span_peaks, key=lambda r: (r.get("sample_id"), r.get("caller"), r.get("peaks_path"))):
        if (row.get("caller") or "").lower() != chosen_caller:
            continue
        sample_id = row.get("sample_id")
        if sample_id and sample_id not in peaks_by_sample:
            peaks_by_sample[sample_id] = row.get("peaks_path")

    missing_peaks = [row.get("sample_id") for row in samples_group if row.get("sample_id") not in peaks_by_sample]
    if missing_peaks:
        reason = "missing_peaks"
        write_failure_outputs(args.outdir, args.group, treated, control, chosen_mode, signature, reason, allow_partial)
        if allow_partial:
            sys.exit(0)
        print(f"ERROR: Missing peaks for samples: {','.join(missing_peaks)}", file=sys.stderr)
        sys.exit(1)

    samples_sorted = sorted(samples_group, key=lambda r: (r.get("condition"), r.get("replicate"), r.get("sample_id")))
    treated_samples = [row for row in samples_sorted if row.get("condition") == treated]
    control_samples = [row for row in samples_sorted if row.get("condition") == control]
    if not treated_samples or not control_samples:
        reason = "missing_condition"
        write_failure_outputs(args.outdir, args.group, treated, control, chosen_mode, signature, reason, allow_partial)
        if allow_partial:
            sys.exit(0)
        print("ERROR: Missing treated/control samples.", file=sys.stderr)
        sys.exit(1)

    spikein_scale_values, spikein_size_factors, use_spikein_factors = compute_spikein_factors(samples_sorted, use_spikein)
    bam_size_factors = {}
    for row, size_factor in zip(samples_sorted, spikein_size_factors):
        bam = row.get("final_bam")
        if bam:
            bam_size_factors[bam] = size_factor if size_factor is not None else 1.0

    diff_path = os.path.join(args.outdir, "span.differential.tsv")
    bed_path = os.path.join(args.outdir, "span.differential.peaks.bed")
    up_path = os.path.join(args.outdir, "span.up.bed")
    down_path = os.path.join(args.outdir, "span.down.bed")
    summary_path = os.path.join(args.outdir, "span.summary.tsv")
    mode_path = os.path.join(args.outdir, "span.mode.txt")
    pooling_path = os.path.join(args.outdir, "span_diff_target_pooling.tsv")

    def finalize_success(rows: List[Dict[str, str]], mode: str, status: str = "RUN", reason: str = "ok"):
        n_tested = len(rows)
        n_fdr = 0
        n_up = 0
        n_down = 0
        for row in rows:
            fdr_val = safe_float(row.get("FDR"))
            log_val = safe_float(row.get("log2FC"))
            if fdr_val is None or log_val is None:
                continue
            if fdr_val <= args.fdr:
                n_fdr += 1
                if log_val > 0:
                    n_up += 1
                elif log_val < 0:
                    n_down += 1
        write_summary(summary_path, args.group, treated, control, n_tested, n_fdr, n_up, n_down, mode, status, reason)

    try:
        if chosen_mode == "native":
            if not native_supported:
                reason = "native_not_supported"
                write_failure_outputs(args.outdir, args.group, treated, control, chosen_mode, signature, reason, allow_partial)
                if allow_partial:
                    sys.exit(0)
                print("ERROR: Native SPAN compare not supported.", file=sys.stderr)
                sys.exit(1)

            flags = parse_flags(help_text)
            treat_flag = flags.get("treatment") or "-t"
            ctrl_flag = flags.get("control") or "-c"
            chrom_flag = flags.get("chrom") or "--cs"
            bin_flag = flags.get("bin") or "--bin"
            gap_flag = flags.get("gap") or "--gap"
            fdr_flag = flags.get("fdr") or "--fdr"
            out_flag = flags.get("output") or "-p"
            out_is_prefix = flags.get("output_is_prefix") == "true"

            if not args.chrom_sizes or not os.path.exists(args.chrom_sizes):
                reason = "missing_chrom_sizes"
                write_failure_outputs(args.outdir, args.group, treated, control, chosen_mode, signature, reason, allow_partial)
                if allow_partial:
                    sys.exit(0)
                print("ERROR: Chrom sizes file required for native SPAN compare.", file=sys.stderr)
                sys.exit(1)

            treated_bams = [row.get("final_bam") for row in treated_samples]
            control_bams = [row.get("final_bam") for row in control_samples]
            for bam in treated_bams + control_bams:
                if not bam or not os.path.exists(bam):
                    raise RuntimeError(f"Missing BAM: {bam}")

            use_lists = signature == "compare:list"
            pooled_rows = []

            def pool_bams(condition: str, bams: List[str]) -> str:
                pooled_name = f"{args.group}.{condition}.pooled.bam"
                pooled_path = os.path.join(args.outdir, pooled_name)
                sorted_bams = sorted(bams)
                cmd = ["samtools", "merge", "-f", "-@", str(max(args.cpus, 1)), pooled_path] + sorted_bams
                run_cmd(cmd, check=True)
                run_cmd(["samtools", "index", pooled_path], check=True)
                if use_spikein_factors:
                    size_vals = [bam_size_factors.get(bam, 1.0) for bam in sorted_bams]
                    pooled_size = sum(size_vals) / len(size_vals) if size_vals else 1.0
                    bam_size_factors[pooled_path] = pooled_size
                pooled_manifest_path = pooled_path
                if args.pooling_dir:
                    pooled_manifest_path = os.path.join(args.pooling_dir, pooled_name)
                pooled_rows.append({
                    "group": args.group,
                    "condition": condition,
                    "pooled_bam": pooled_manifest_path,
                    "source_bams": ",".join(sorted_bams),
                })
                return pooled_path

            treated_arg = ""
            control_arg = ""
            compare_treated_bams = []
            compare_control_bams = []
            if use_lists:
                treated_arg = ",".join(treated_bams)
                control_arg = ",".join(control_bams)
                compare_treated_bams = treated_bams
                compare_control_bams = control_bams
            else:
                if len(treated_bams) > 1:
                    treated_arg = pool_bams(treated, treated_bams)
                else:
                    treated_arg = treated_bams[0]
                if len(control_bams) > 1:
                    control_arg = pool_bams(control, control_bams)
                else:
                    control_arg = control_bams[0]
                compare_treated_bams = [treated_arg]
                compare_control_bams = [control_arg]

            output_target = "span.native"
            output_path = output_target
            if out_is_prefix:
                output_path = output_target
            else:
                output_path = f"{output_target}.peak"

            cmd = [
                "java",
                f"-Xmx{args.java_heap}",
                "-jar",
                args.jar,
                "compare",
                treat_flag,
                treated_arg,
                ctrl_flag,
                control_arg,
                chrom_flag,
                args.chrom_sizes,
                bin_flag,
                str(args.bin),
                gap_flag,
                str(args.gap),
                fdr_flag,
                str(args.fdr),
                out_flag,
                output_path,
            ]

            result = run_cmd(cmd, check=False)
            if result.returncode != 0:
                raise RuntimeError(result.stderr.strip() or "native_compare_failed")

            native_output = None
            if out_is_prefix:
                candidates = [f for f in os.listdir(args.outdir) if f.startswith(output_target)]
                candidates = [c for c in candidates if c.endswith(".peak") or c.endswith(".bed") or c.endswith(".tsv")]
                candidates.sort()
                if candidates:
                    native_output = os.path.join(args.outdir, candidates[0])
            else:
                native_output = os.path.join(args.outdir, output_path)

            if not native_output or not os.path.exists(native_output):
                raise RuntimeError("native_output_missing")

            rows = parse_native_output(native_output)
            if not rows:
                raise RuntimeError("native_output_empty")

            bed_rows = [(r["chr"], int(r["start"]), int(r["end"])) for r in rows if is_int(r["start"]) and is_int(r["end"])]
            bed_rows.sort(key=lambda x: (x[0], x[1], x[2]))
            write_bed(bed_path, bed_rows)

            pseudocount = 0.5
            treated_size_factors = None
            control_size_factors = None
            if use_spikein_factors:
                treated_size_factors = [bam_size_factors.get(bam, 1.0) for bam in compare_treated_bams]
                control_size_factors = [bam_size_factors.get(bam, 1.0) for bam in compare_control_bams]
            orientation_map = compute_orientation_log2fc(
                bed_path,
                compare_treated_bams,
                compare_control_bams,
                pseudocount,
                treated_size_factors,
                control_size_factors,
            )
            for row in rows:
                key = (row.get("chr"), row.get("start"), row.get("end"))
                values = orientation_map.get(key, [])
                log2fc_user = values.pop(0) if values else None
                row["log2FC_span"] = row.get("log2FC", "NA")
                row["log2FC"] = f"{log2fc_user:.6f}" if log2fc_user is not None else "NA"

            write_tsv(diff_path, ["chr", "start", "end", "log2FC", "log2FC_span", "pval", "FDR"], rows)

            up_rows = filter_bed_by_stats(rows, args.fdr, "up")
            down_rows = filter_bed_by_stats(rows, args.fdr, "down")
            up_rows.sort(key=lambda x: (x[0], x[1], x[2]))
            down_rows.sort(key=lambda x: (x[0], x[1], x[2]))
            write_bed(up_path, up_rows)
            write_bed(down_path, down_rows)

            normalization_path = os.path.join(args.outdir, "span.normalization_factors.tsv")
            write_normalization_factors(
                normalization_path,
                samples_sorted,
                spikein_scale_values,
                spikein_size_factors,
                use_spikein_factors,
                "native",
            )
            print(f"INFO: wrote normalization factors to {normalization_path}", file=sys.stderr)

            if pooled_rows:
                write_tsv(
                    pooling_path,
                    ["group", "condition", "pooled_bam", "source_bams"],
                    pooled_rows,
                )

            write_mode(
                mode_path,
                "native",
                signature,
                notes=[f"orientation=treated/control", f"orientation_pseudocount={pseudocount}"],
            )
            finalize_success(rows, "native")
            return

        # Fallback path
        all_intervals = []
        for row in samples_sorted:
            peak_path = peaks_by_sample.get(row.get("sample_id"))
            if not peak_path or not os.path.exists(peak_path):
                raise RuntimeError(f"Missing peaks file: {peak_path}")
            all_intervals.extend(read_intervals(peak_path))

        merged = merge_intervals(all_intervals)
        if not merged:
            raise RuntimeError("no_regions")

        write_bed(bed_path, merged)

        # counts
        sample_ids = [row.get("sample_id") for row in samples_sorted]
        bams = [row.get("final_bam") for row in samples_sorted]
        for bam in bams:
            if not bam or not os.path.exists(bam):
                raise RuntimeError(f"Missing BAM: {bam}")

        counts_path = os.path.join(args.outdir, "span.counts.tsv")
        cmd = ["bedtools", "multicov", "-bed", bed_path, "-bams"] + bams
        result = run_cmd(cmd, check=False)
        if result.returncode != 0:
            raise RuntimeError(result.stderr.strip() or "bedtools_multicov_failed")
        with open(counts_path, "w") as handle:
            handle.write("chr\tstart\tend\t" + "\t".join(sample_ids) + "\n")
            handle.write(result.stdout)

        samples_meta_path = os.path.join(args.outdir, "span.samples.tsv")
        with open(samples_meta_path, "w") as handle:
            handle.write("sample_id\tcondition\tspikein_scale_factor\n")
            for row in samples_sorted:
                handle.write(
                    f"{row.get('sample_id','')}\t{row.get('condition','')}\t{row.get('spikein_scale_factor','NA')}\n"
                )

        normalization_path = os.path.join(args.outdir, "span.normalization_factors.tsv")
        write_normalization_factors(
            normalization_path,
            samples_sorted,
            spikein_scale_values,
            spikein_size_factors,
            use_spikein_factors,
            args.fallback_backend,
        )
        print(f"INFO: wrote normalization factors to {normalization_path}", file=sys.stderr)

        r_script = os.path.join(os.path.dirname(__file__), "span_fallback_de.R")
        cmd = [
            "Rscript",
            r_script,
            "--counts",
            counts_path,
            "--samples",
            samples_meta_path,
            "--contrast",
            args.contrast,
            "--backend",
            args.fallback_backend,
            "--use-spikein",
            "true" if use_spikein else "false",
            "--out",
            diff_path,
        ]
        result = run_cmd(cmd, check=False)
        if result.returncode != 0:
            raise RuntimeError(result.stderr.strip() or "fallback_de_failed")

        rows = parse_native_output(diff_path)
        if not rows:
            raise RuntimeError("fallback_output_empty")

        up_rows = filter_bed_by_stats(rows, args.fdr, "up")
        down_rows = filter_bed_by_stats(rows, args.fdr, "down")
        up_rows.sort(key=lambda x: (x[0], x[1], x[2]))
        down_rows.sort(key=lambda x: (x[0], x[1], x[2]))
        write_bed(up_path, up_rows)
        write_bed(down_path, down_rows)

        write_mode(mode_path, "fallback", signature)
        finalize_success(rows, "fallback")
    except Exception as exc:
        reason = str(exc).split("\n")[0] if exc else "failed"
        write_failure_outputs(args.outdir, args.group, treated, control, chosen_mode, signature, reason, allow_partial)
        if allow_partial:
            sys.exit(0)
        print(f"ERROR: SPAN diff failed: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
