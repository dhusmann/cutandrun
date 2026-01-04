#!/usr/bin/env python3
import argparse
import csv
import os
import subprocess
import sys


def read_tsv(path):
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for row in reader]


def detect_native_support(jar_path, java_heap):
    try:
        result = subprocess.run(
            ["java", f"-Xmx{java_heap}", "-jar", jar_path, "--help"],
            capture_output=True,
            text=True,
            check=False,
        )
        output = (result.stdout or "") + (result.stderr or "")
    except FileNotFoundError:
        return False, "java_missing"
    except Exception:
        return False, "help_failed"
    if "compare" in output.lower():
        return True, "compare"
    return False, "no_compare"


def read_peaks(peaks_path, max_rows=200):
    regions = []
    with open(peaks_path, "r") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            regions.append((parts[0], parts[1], parts[2]))
            if len(regions) >= max_rows:
                break
    return regions


def main():
    parser = argparse.ArgumentParser(description="Stub SPAN/OmniPeak differential")
    parser.add_argument("--jar", required=True)
    parser.add_argument("--mode", required=True)
    parser.add_argument("--contrast", required=True)
    parser.add_argument("--group", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--peaks", required=True)
    parser.add_argument("--bin", type=int, default=200)
    parser.add_argument("--gap", type=int, default=5)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--fallback-backend", default="DESeq2")
    parser.add_argument("--java-heap", default="8G")
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()

    contrast = [x.strip() for x in args.contrast.split(",") if x.strip()]
    treated = contrast[0] if len(contrast) > 0 else "treated"
    control = contrast[1] if len(contrast) > 1 else "control"

    native_supported, signature = detect_native_support(args.jar, args.java_heap)
    chosen_mode = "fallback"
    status = "RUN"
    reason = "ok"
    if args.mode == "native":
        if not native_supported:
            if signature == "java_missing":
                print("WARNING: Java not found; falling back to SPAN fallback mode.", file=sys.stderr)
                chosen_mode = "fallback"
                reason = "java_missing_fallback"
            else:
                print("ERROR: SPAN jar does not support native compare.", file=sys.stderr)
                sys.exit(1)
        else:
            chosen_mode = "native"
    elif args.mode == "fallback":
        chosen_mode = "fallback"
    else:
        chosen_mode = "native" if native_supported else "fallback"
        if not native_supported and signature == "java_missing":
            reason = "java_missing_fallback"

    if not os.path.exists(args.peaks):
        print("ERROR: SPAN peaks file not found for fallback.", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    regions = read_peaks(args.peaks)

    diff_path = os.path.join(args.outdir, "span.differential.tsv")
    bed_path = os.path.join(args.outdir, "span.differential.peaks.bed")
    up_path = os.path.join(args.outdir, "span.up.bed")
    down_path = os.path.join(args.outdir, "span.down.bed")
    summary_path = os.path.join(args.outdir, "span.summary.tsv")
    mode_path = os.path.join(args.outdir, "span.mode.txt")

    with open(diff_path, "w") as handle:
        handle.write("chr\tstart\tend\tlog2FC\tpval\tFDR\n")
        for chrom, start, end in regions:
            handle.write(f"{chrom}\t{start}\t{end}\t0.0\t1.0\t1.0\n")

    with open(bed_path, "w") as handle:
        for chrom, start, end in regions:
            handle.write(f"{chrom}\t{start}\t{end}\n")

    open(up_path, "w").close()
    open(down_path, "w").close()

    with open(summary_path, "w") as handle:
        handle.write("group\ttreated\tcontrol\tn_tested\tn_fdr_pass\tn_up\tn_down\tmode\tstatus\treason\n")
        handle.write(f"{args.group}\t{treated}\t{control}\t{len(regions)}\t0\t0\t0\t{chosen_mode}\t{status}\t{reason}\n")

    with open(mode_path, "w") as handle:
        handle.write(f"{chosen_mode}:{signature}\n")


if __name__ == "__main__":
    main()
