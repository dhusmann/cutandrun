#!/usr/bin/env python3
import argparse
import csv
import os
import sys


def detect_delimiter(path):
    with open(path, "r", newline="") as handle:
        header = handle.readline()
    if "\t" in header:
        return "\t"
    return ","


def read_table(path):
    delim = detect_delimiter(path)
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delim)
        rows = [row for row in reader]
    return rows


def load_ms_coeffs(path):
    if not path:
        return {}
    rows = read_table(path)
    if not rows:
        return {}
    header = rows[0].keys()
    if "sample_id" not in header or "ms_coeff" not in header:
        raise ValueError("ms_coeffs file must contain columns: sample_id, ms_coeff")
    coeffs = {}
    for row in rows:
        sample_id = row.get("sample_id")
        if not sample_id:
            continue
        coeffs[sample_id] = row.get("ms_coeff") or "NA"
    return coeffs


def sort_key(value):
    if value is None:
        return ""
    return str(value)


def write_samples(samples_rows, ms_coeffs, normalisation_mode, out_path):
    header = [
        "sample_id",
        "group",
        "condition",
        "replicate",
        "final_bam",
        "final_bai",
        "normalisation_mode",
        "spikein_scale_factor",
        "ms_coeff",
        "bigwig_path",
        "input_bam",
        "input_bai",
    ]
    rows = []
    for row in samples_rows:
        sample_id = row.get("sample_id") or ""
        group = row.get("group") or ""
        condition = row.get("condition") or ""
        replicate = row.get("replicate") or ""
        final_bam = row.get("final_bam") or "NA"
        final_bai = row.get("final_bai") or "NA"
        spikein = row.get("spikein_scale_factor") or "NA"
        bigwig = row.get("bigwig_path") or "NA"
        input_bam = row.get("input_bam") or "NA"
        input_bai = row.get("input_bai") or "NA"
        ms_coeff = ms_coeffs.get(sample_id, "NA")
        rows.append({
            "sample_id": sample_id,
            "group": group,
            "condition": condition,
            "replicate": replicate,
            "final_bam": final_bam,
            "final_bai": final_bai,
            "normalisation_mode": normalisation_mode,
            "spikein_scale_factor": spikein,
            "ms_coeff": ms_coeff,
            "bigwig_path": bigwig,
            "input_bam": input_bam,
            "input_bai": input_bai,
        })

    rows.sort(key=lambda r: (r["group"], "NA", r["condition"], sort_key(r["replicate"]), r["sample_id"]))

    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_peaks(peaks_rows, out_path):
    header = [
        "sample_id",
        "group",
        "condition",
        "caller",
        "peaks_path",
        "peaks_format",
    ]
    rows = []
    for row in peaks_rows:
        rows.append({
            "sample_id": row.get("sample_id") or "",
            "group": row.get("group") or "",
            "condition": row.get("condition") or "",
            "caller": row.get("caller") or "",
            "peaks_path": row.get("peaks_path") or "NA",
            "peaks_format": row.get("peaks_format") or "NA",
            "replicate": row.get("replicate") or "",
        })

    rows.sort(key=lambda r: (r["group"], r["caller"], r["condition"], sort_key(r["replicate"]), r["sample_id"]))

    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description="Build differential manifests")
    parser.add_argument("--samples", required=True)
    parser.add_argument("--peaks", required=True)
    parser.add_argument("--normalisation-mode", required=True)
    parser.add_argument("--ms-coeffs")
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()

    samples_rows = read_table(args.samples)
    peaks_rows = read_table(args.peaks)
    ms_coeffs = {}
    if args.ms_coeffs:
        try:
            ms_coeffs = load_ms_coeffs(args.ms_coeffs)
        except Exception as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            sys.exit(1)

    out_samples = os.path.join(args.outdir, "differential_manifest.samples.tsv")
    out_peaks = os.path.join(args.outdir, "differential_manifest.peaks.tsv")

    write_samples(samples_rows, ms_coeffs, args.normalisation_mode, out_samples)
    write_peaks(peaks_rows, out_peaks)


if __name__ == "__main__":
    main()
