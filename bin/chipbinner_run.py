#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import random


def read_tsv(path):
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for row in reader]


def read_chrom_sizes(path):
    chroms = []
    with open(path, "r") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            chroms.append((parts[0], int(parts[1])))
    return chroms


def write_windows(chroms, bin_size, out_path, max_bins=200):
    bins = []
    for chrom, size in chroms:
        start = 0
        while start < size and len(bins) < max_bins:
            end = min(start + bin_size, size)
            bins.append((chrom, start, end))
            start = end
        if len(bins) >= max_bins:
            break
    with open(out_path, "w") as handle:
        for chrom, start, end in bins:
            handle.write(f"{chrom}\t{start}\t{end}\n")
    return bins


def main():
    parser = argparse.ArgumentParser(description="Stub ChIPBinner runner")
    parser.add_argument("--samples", required=True)
    parser.add_argument("--group", required=True)
    parser.add_argument("--contrast", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--bin-size", type=int, required=True)
    parser.add_argument("--windows-dir")
    parser.add_argument("--use-input", action="store_true")
    parser.add_argument("--pseudocount", type=int, default=1)
    parser.add_argument("--grid-minpts", default="")
    parser.add_argument("--grid-minsamps", default="")
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--lfc", type=float, default=1.0)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--k-value", type=int, default=100000)
    parser.add_argument("--functional-db")
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    plots_dir = os.path.join(args.outdir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    samples = read_tsv(args.samples)
    contrast = [x.strip() for x in args.contrast.split(",") if x.strip()]
    treated = contrast[0] if len(contrast) > 0 else "treated"
    control = contrast[1] if len(contrast) > 1 else "control"

    # Write samplesheet
    samplesheet_path = os.path.join(args.outdir, "chipbinner.samplesheet.csv")
    with open(samplesheet_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", "group", "condition", "replicate", "bam", "bigwig", "spikein_scale_factor", "ms_coeff"])
        for row in samples:
            writer.writerow([
                row.get("sample_id", ""),
                row.get("group", ""),
                row.get("condition", ""),
                row.get("replicate", ""),
                row.get("final_bam", ""),
                row.get("bigwig_path", ""),
                row.get("spikein_scale_factor", ""),
                row.get("ms_coeff", ""),
            ])

    # Windows
    chroms = read_chrom_sizes(args.chrom_sizes)
    windows_path = os.path.join(args.outdir, "chipbinner.windows.bed")
    bins = []
    if args.windows_dir:
        if os.path.isfile(args.windows_dir):
            windows_path = args.windows_dir
        elif os.path.isdir(args.windows_dir):
            candidates = sorted(glob.glob(os.path.join(args.windows_dir, "*.bed")))
            if candidates:
                windows_path = candidates[0]
    if os.path.exists(windows_path) and os.path.abspath(windows_path) != os.path.abspath(os.path.join(args.outdir, "chipbinner.windows.bed")):
        with open(windows_path, "r") as handle:
            for line in handle:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    bins.append((parts[0], int(parts[1]), int(parts[2])))
        with open(os.path.join(args.outdir, "chipbinner.windows.bed"), "w") as out_handle:
            for chrom, start, end in bins:
                out_handle.write(f"{chrom}\t{start}\t{end}\n")
        windows_path = os.path.join(args.outdir, "chipbinner.windows.bed")
    if not bins:
        bins = write_windows(chroms, args.bin_size, windows_path)

    # Counts and normalized matrix
    count_path = os.path.join(args.outdir, "chipbinner.bin_counts.tsv")
    norm_path = os.path.join(args.outdir, "chipbinner.normalized_matrix.tsv")
    sample_ids = [row.get("sample_id", "") for row in samples]
    with open(count_path, "w") as handle:
        handle.write("bin\t" + "\t".join(sample_ids) + "\n")
        for idx, (chrom, start, end) in enumerate(bins):
            values = [str((idx + 1) * (s_idx + 1)) for s_idx in range(len(sample_ids))]
            handle.write(f"{chrom}:{start}-{end}\t" + "\t".join(values) + "\n")
    with open(norm_path, "w") as handle:
        handle.write("bin\t" + "\t".join(sample_ids) + "\n")
        for idx, (chrom, start, end) in enumerate(bins):
            values = [str((idx + 1) * (s_idx + 1)) for s_idx in range(len(sample_ids))]
            handle.write(f"{chrom}:{start}-{end}\t" + "\t".join(values) + "\n")

    # HDBSCAN grid summary
    grid_path = os.path.join(args.outdir, "chipbinner.hdbscan_grid_summary.tsv")
    minpts_values = [x for x in args.grid_minpts.split(",") if x] or ["100"]
    minsamps_values = [x for x in args.grid_minsamps.split(",") if x] or ["100"]
    with open(grid_path, "w") as handle:
        handle.write("minPts\tminSamps\tn_clusters\tstatus\n")
        for mpt in minpts_values:
            for ms in minsamps_values:
                handle.write(f"{mpt}\t{ms}\t2\tOK\n")

    # Clusters
    clusters_path = os.path.join(args.outdir, "chipbinner.clusters.tsv")
    with open(clusters_path, "w") as handle:
        handle.write("chrom\tstart\tend\tcluster\n")
        for idx, (chrom, start, end) in enumerate(bins):
            handle.write(f"{chrom}\t{start}\t{end}\t1\n")

    # Differential
    diff_path = os.path.join(args.outdir, "chipbinner.differential.tsv")
    with open(diff_path, "w") as handle:
        handle.write("chr\tstart\tend\tlog2FC\tpval\tFDR\tcluster\n")
        for idx, (chrom, start, end) in enumerate(bins):
            handle.write(f"{chrom}\t{start}\t{end}\t0.0\t1.0\t1.0\t1\n")

    # Summary
    summary_path = os.path.join(args.outdir, "chipbinner.summary.tsv")
    with open(summary_path, "w") as handle:
        handle.write("group\ttreated\tcontrol\tn_bins_tested\tn_fdr_pass\tn_up\tn_down\tn_clusters\tchosen_minPts\tchosen_minSamps\tstatus\treason\n")
        handle.write(f"{args.group}\t{treated}\t{control}\t{len(bins)}\t0\t0\t0\t1\t{minpts_values[0]}\t{minsamps_values[0]}\tRUN\tok\n")

    # Enrichment stub
    enrichment_dir = os.path.join(args.outdir, "enrichment")
    os.makedirs(enrichment_dir, exist_ok=True)
    enrichment_file = os.path.join(enrichment_dir, "not_run.txt")
    if args.functional_db:
        enrichment_file = os.path.join(enrichment_dir, "enrichment_stub.txt")
    with open(enrichment_file, "w") as handle:
        handle.write("Enrichment not run\n" if not args.functional_db else "Enrichment stub\n")

    with open(os.path.join(plots_dir, "chipbinner_stub.txt"), "w") as handle:
        handle.write("stub plots\n")


if __name__ == "__main__":
    main()
