#!/usr/bin/env python3
import argparse
import csv
import glob
import hashlib
import math
import os
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

try:
    import hdbscan
except ImportError as exc:
    hdbscan = None

try:
    from scipy.stats import fisher_exact
except ImportError:
    fisher_exact = None


def parse_args():
    parser = argparse.ArgumentParser(description="Run ChIPBinner workflow")
    parser.add_argument("--samples", required=True)
    parser.add_argument("--group", required=True)
    parser.add_argument("--contrast", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--bin-size", type=int, required=True)
    parser.add_argument("--windows")
    parser.add_argument("--windows-dir")
    parser.add_argument("--blacklist")
    parser.add_argument("--use-input", action="store_true")
    parser.add_argument("--use-spikein", action="store_true")
    parser.add_argument("--pseudocount", type=float, default=1.0)
    parser.add_argument("--grid-minpts", default="")
    parser.add_argument("--grid-minsamps", default="")
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--lfc", type=float, default=1.0)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--k-value", type=int, default=100000)
    parser.add_argument("--functional-db")
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument("--outdir", default=".")
    return parser.parse_args()


def read_tsv(path):
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for row in reader]


def parse_float(value, default=None):
    if value is None:
        return default
    text = str(value).strip()
    if text in ("", "NA", "None", "none"):
        return default
    try:
        return float(text)
    except Exception:
        return default


def genome_id_from_path(path):
    base = os.path.basename(path)
    if base.endswith(".sizes"):
        base = base[: -len(".sizes")]
    else:
        base = os.path.splitext(base)[0]
    return base or "genome"


def hash_file(path):
    if not path or not os.path.exists(path):
        return None
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()[:12]


def ensure_dir(path):
    if path:
        os.makedirs(path, exist_ok=True)


def run_cmd(cmd, stdout=None):
    result = subprocess.run(cmd, stdout=stdout, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({' '.join(cmd)}): {result.stderr.strip()}")


def choose_windows_file(windows_path, windows_dir, bin_size, chrom_sizes=None, blacklist=None):
    if windows_path:
        return windows_path
    if not windows_dir:
        return None
    if os.path.isfile(windows_dir):
        return windows_dir
    if os.path.isdir(windows_dir):
        genome_id = genome_id_from_path(chrom_sizes) if chrom_sizes else None
        blacklist_hash = hash_file(blacklist) if blacklist else None
        preferred = []
        if genome_id:
            if blacklist_hash:
                preferred.append(f"windows.{genome_id}.{bin_size}.{blacklist_hash}.bed")
            preferred.append(f"windows.{genome_id}.{bin_size}.bed")
        else:
            preferred.append(f"windows.{bin_size}.bed")
        for name in preferred:
            candidate = os.path.join(windows_dir, name)
            if os.path.exists(candidate):
                return candidate
        raise RuntimeError(
            f"No windows BED found in {windows_dir}. Expected one of: {', '.join(preferred)}"
        )
    return None


def generate_windows(chrom_sizes, bin_size, blacklist, out_path):
    with tempfile.TemporaryDirectory() as tmpdir:
        raw_path = os.path.join(tmpdir, "windows_raw.bed")
        with open(raw_path, "w") as handle:
            with open(chrom_sizes, "r") as chrom_handle:
                for line in chrom_handle:
                    parts = line.strip().split("\t")
                    if len(parts) < 2:
                        continue
                    chrom, size = parts[0], int(parts[1])
                    start = 0
                    while start < size:
                        end = min(start + bin_size, size)
                        handle.write(f"{chrom}\t{start}\t{end}\n")
                        start = end
        work_path = raw_path
        if blacklist and os.path.exists(blacklist):
            filtered_path = os.path.join(tmpdir, "windows_filtered.bed")
            run_cmd(["bedtools", "subtract", "-a", raw_path, "-b", blacklist], stdout=open(filtered_path, "w"))
            work_path = filtered_path
        with open(out_path, "w") as out_handle:
            run_cmd(["sort", "-k1,1", "-k2,2n", work_path], stdout=out_handle)


def load_windows(path):
    windows = pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end"])
    windows["bin_id"] = windows.apply(lambda r: f"{r.chrom}:{int(r.start)}-{int(r.end)}", axis=1)
    return windows


def write_samplesheet(samples, out_path):
    header = [
        "sample_id",
        "group",
        "condition",
        "replicate",
        "final_bam",
        "final_bai",
        "bigwig_path",
        "spikein_scale_factor",
        "ms_coeff",
        "input_bam",
        "input_bai",
    ]
    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        writer.writeheader()
        for row in samples:
            writer.writerow({
                "sample_id": row.get("sample_id", ""),
                "group": row.get("group", ""),
                "condition": row.get("condition", ""),
                "replicate": row.get("replicate", ""),
                "final_bam": row.get("final_bam", ""),
                "final_bai": row.get("final_bai", ""),
                "bigwig_path": row.get("bigwig_path", ""),
                "spikein_scale_factor": row.get("spikein_scale_factor", ""),
                "ms_coeff": row.get("ms_coeff", ""),
                "input_bam": row.get("input_bam", ""),
                "input_bai": row.get("input_bai", ""),
            })

def compute_counts(windows_path, samples, out_path, use_input=False):
    bam_paths = []
    sample_ids = []
    for row in samples:
        bam = row.get("final_bam")
        if not bam or bam in ("NA", "None"):
            raise RuntimeError(f"Missing final_bam for sample {row.get('sample_id')}")
        if not os.path.exists(bam):
            raise RuntimeError(f"BAM not found: {bam}")
        bam_paths.append(bam)
        sample_ids.append(row.get("sample_id", os.path.basename(bam)))

    cmd = ["bedtools", "multicov", "-bed", windows_path, "-bams"] + bam_paths
    with open(out_path, "w") as handle:
        run_cmd(cmd, stdout=handle)

    cols = ["chrom", "start", "end"] + sample_ids
    counts = pd.read_csv(out_path, sep="\t", header=None, names=cols)
    counts["bin_id"] = counts.apply(lambda r: f"{r.chrom}:{int(r.start)}-{int(r.end)}", axis=1)

    if use_input:
        input_counts_map = {}
        for row in samples:
            bam = row.get("input_bam")
            if not bam or bam in ("NA", "None", "none", ""):
                raise RuntimeError("chipbinner_use_input requested but no input_bam values were provided")
            if bam in input_counts_map:
                continue
            if not os.path.exists(bam):
                raise RuntimeError(f"Input BAM not found: {bam}")
            with tempfile.NamedTemporaryFile(delete=False) as tmp_handle:
                tmp_path = tmp_handle.name
            cmd = ["bedtools", "multicov", "-bed", windows_path, "-bams", bam]
            with open(tmp_path, "w") as handle:
                run_cmd(cmd, stdout=handle)
            tmp_df = pd.read_csv(tmp_path, sep="\t", header=None, names=["chrom", "start", "end", "count"])
            os.unlink(tmp_path)
            input_counts_map[bam] = tmp_df["count"].astype(float).values

        for sample_id, row in zip(sample_ids, samples):
            bam = row.get("input_bam")
            counts[sample_id] = counts[sample_id].astype(float).sub(input_counts_map[bam], axis=0).clip(lower=0)
    return counts, sample_ids


def write_matrix(counts_df, sample_ids, out_path):
    matrix = counts_df[["chrom", "start", "end", "bin_id"] + sample_ids]
    matrix.to_csv(out_path, sep="\t", index=False)


def normalize_counts(counts_df, sample_ids, samples, use_spikein, pseudocount):
    counts = counts_df[sample_ids].astype(float)
    if use_spikein:
        missing = []
        for row in samples:
            spike_raw = row.get("spikein_scale_factor")
            if parse_float(spike_raw) is None:
                missing.append(row.get("sample_id", ""))
        if missing:
            raise RuntimeError(f"Missing spike-in scale factors for samples: {', '.join(sorted(set(missing)))}")
    factors = []
    scaling_applied = []
    for idx, row in enumerate(samples):
        steps = []
        applied = False
        spike_raw = row.get("spikein_scale_factor")
        spike = parse_float(spike_raw) if use_spikein else None
        spike_size = None
        if spike is not None:
            if spike == 0:
                raise RuntimeError("Spike-in scale factor cannot be zero")
            spike_size = 1.0 / spike
            # Apply the scale factor directly (equivalent to dividing by size factor).
            counts.iloc[:, idx] = counts.iloc[:, idx] * spike
            steps.append("spikein")
            applied = True

        ms_raw = row.get("ms_coeff")
        ms_coeff = parse_float(ms_raw)
        ms_size = None
        if ms_coeff is not None:
            if ms_coeff == 0:
                raise RuntimeError("MS coefficient cannot be zero")
            ms_size = 1.0 / ms_coeff
            counts.iloc[:, idx] = counts.iloc[:, idx] / ms_size
            steps.append("ms_coeff")
            applied = True

        factors.append({
            "sample_id": row.get("sample_id", ""),
            "group": row.get("group", ""),
            "condition": row.get("condition", ""),
            "spikein_scale_factor": spike if spike is not None else "NA",
            "spikein_size_factor": spike_size if spike_size is not None else "NA",
            "ms_coeff": ms_coeff if ms_coeff is not None else "NA",
            "ms_size_factor": ms_size if ms_size is not None else "NA",
            "applied_steps": steps,
        })
        scaling_applied.append(applied)

    apply_cpm = not any(scaling_applied)
    if apply_cpm:
        library_size = counts.sum(axis=0)
        library_size[library_size == 0] = 1.0
        normalized = counts.divide(library_size, axis=1) * 1e6
    else:
        normalized = counts.copy()
    normalized = normalized + pseudocount

    for entry in factors:
        steps = entry.get("applied_steps", [])
        steps.append("pseudocount")
        if apply_cpm:
            steps.append("library_size_cpm")
        entry["applied_steps"] = ",".join(steps)

    return normalized, pd.DataFrame(factors)


def plot_pca(norm_matrix, sample_ids, samples, out_path):
    data = np.log2(norm_matrix + 1)
    data = data.T
    labels = [row.get("condition", "") for row in samples]
    pca = PCA(n_components=2)
    coords = pca.fit_transform(data)
    plt.figure(figsize=(6, 5))
    for condition in sorted(set(labels)):
        idxs = [i for i, c in enumerate(labels) if c == condition]
        plt.scatter(coords[idxs, 0], coords[idxs, 1], label=condition)
    for i, sample_id in enumerate(sample_ids):
        plt.text(coords[i, 0], coords[i, 1], sample_id, fontsize=7)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend(loc="best", fontsize=7)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def plot_correlation(norm_matrix, sample_ids, out_path):
    data = np.log2(norm_matrix + 1)
    corr = pd.DataFrame(data, columns=sample_ids).corr()
    plt.figure(figsize=(6, 5))
    sns.heatmap(corr, cmap="vlag", center=0, xticklabels=True, yticklabels=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def parse_grid(values, default_values):
    if values:
        items = [v.strip() for v in str(values).split(",") if v.strip()]
        return [int(x) for x in items]
    return default_values


def hdbscan_grid(norm_matrix, windows_df, minpts_list, minsamps_list, out_dir, summary_path):
    if hdbscan is None:
        raise RuntimeError("hdbscan python package is required for clustering")
    ensure_dir(out_dir)
    features = StandardScaler().fit_transform(norm_matrix.values)
    results = []
    assignments = {}
    n_bins = features.shape[0]
    for minpts in minpts_list:
        for minsamps in minsamps_list:
            status = "OK"
            n_clusters = 0
            noise_fraction = 1.0
            silhouette = np.nan
            labels = None
            if minpts > n_bins:
                status = "INVALID"
            else:
                try:
                    clusterer = hdbscan.HDBSCAN(min_cluster_size=minpts, min_samples=minsamps)
                    labels = clusterer.fit_predict(features)
                    n_clusters = len({lab for lab in labels if lab >= 0})
                    noise_fraction = float(np.sum(labels < 0)) / float(len(labels))
                    if n_clusters > 1:
                        mask = labels >= 0
                        silhouette = silhouette_score(features[mask], labels[mask])
                except Exception as exc:
                    status = f"FAIL:{exc}"
            results.append({
                "minPts": minpts,
                "minSamps": minsamps,
                "n_clusters": n_clusters,
                "noise_fraction": noise_fraction,
                "silhouette": silhouette,
                "status": status,
                "labels": labels,
            })
    # selection rule: maximize silhouette (if available), then minimize noise, then max clusters, then smallest params
    best = None
    for entry in results:
        if entry["status"] != "OK":
            continue
        sil = entry["silhouette"]
        if math.isnan(sil):
            score = -1.0 - entry["noise_fraction"]
        else:
            score = sil - entry["noise_fraction"]
        entry["score"] = score
        if best is None:
            best = entry
            continue
        if score > best["score"]:
            best = entry
            continue
        if score == best["score"]:
            if entry["n_clusters"] > best["n_clusters"]:
                best = entry
                continue
            if entry["n_clusters"] == best["n_clusters"] and entry["noise_fraction"] < best["noise_fraction"]:
                best = entry
                continue
            if entry["n_clusters"] == best["n_clusters"] and entry["noise_fraction"] == best["noise_fraction"]:
                if (entry["minPts"], entry["minSamps"]) < (best["minPts"], best["minSamps"]):
                    best = entry

    if best is None:
        raise RuntimeError("No successful HDBSCAN grid points")

    for entry in results:
        if entry["status"] != "OK" or entry["labels"] is None:
            continue
        labels = entry["labels"]
        out_path = os.path.join(out_dir, f"clusters_minPts{entry['minPts']}_minSamps{entry['minSamps']}.tsv")
        df = windows_df.copy()
        df["cluster"] = labels
        df.to_csv(out_path, sep="\t", index=False)
        assignments[(entry["minPts"], entry["minSamps"])] = out_path

    summary_rows = []
    for entry in results:
        sil = entry["silhouette"]
        sil_val = "NA" if math.isnan(sil) else f"{sil:.4f}"
        score = entry.get("score")
        score_val = "NA" if score is None else f"{score:.4f}"
        selected = "true" if best and entry["minPts"] == best["minPts"] and entry["minSamps"] == best["minSamps"] else "false"
        summary_rows.append({
            "minPts": entry["minPts"],
            "minSamps": entry["minSamps"],
            "n_clusters": entry["n_clusters"],
            "noise_fraction": f"{entry['noise_fraction']:.4f}",
            "silhouette": sil_val,
            "score": score_val,
            "status": entry["status"],
            "selected": selected,
        })

    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)

    chosen = {
        "minPts": best["minPts"],
        "minSamps": best["minSamps"],
        "n_clusters": best["n_clusters"],
        "labels": best["labels"],
        "summary_path": summary_path,
    }
    return chosen, summary_path


def run_rots(matrix_path, samplesheet_path, treated, control, bootstrap, k_value, out_path):
    script_path = os.path.join(os.path.dirname(__file__), "chipbinner_rots.R")
    cmd = [
        "Rscript",
        script_path,
        "--matrix", matrix_path,
        "--samplesheet", samplesheet_path,
        "--treated", treated,
        "--control", control,
        "--bootstrap", str(bootstrap),
        "--k-value", str(k_value),
        "--out", out_path,
    ]
    run_cmd(cmd)


def bh_fdr(pvals):
    pvals = np.asarray(pvals)
    n = len(pvals)
    order = np.argsort(pvals)
    ranks = np.arange(1, n + 1)
    fdr = np.empty(n, dtype=float)
    fdr[order] = pvals[order] * n / ranks
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    fdr[fdr > 1] = 1
    return fdr


def run_enrichment(windows_df, clusters_path, functional_db, out_dir):
    ensure_dir(out_dir)
    if not functional_db:
        stub = os.path.join(out_dir, "not_run.txt")
        with open(stub, "w") as handle:
            handle.write("NOT_RUN\n")
        return "NOT_RUN"
    if fisher_exact is None:
        raise RuntimeError("scipy is required for enrichment")
    if os.path.isdir(functional_db):
        bed_files = sorted(glob.glob(os.path.join(functional_db, "*.bed")))
    else:
        bed_files = [functional_db]
    if not bed_files:
        raise RuntimeError("No enrichment BED files found")

    # write background windows bed
    bg_bed = os.path.join(out_dir, "windows.bed")
    windows_df[["chrom", "start", "end"]].to_csv(bg_bed, sep="\t", header=False, index=False)

    clusters = pd.read_csv(clusters_path, sep="\t")
    clusters = clusters[clusters["cluster"] >= 0]
    if clusters.empty:
        stub = os.path.join(out_dir, "no_clusters.txt")
        with open(stub, "w") as handle:
            handle.write("NO_CLUSTERS\n")
        return "NOT_RUN"

    all_results = []
    bg_total = len(windows_df)
    for cluster_id in sorted(clusters["cluster"].unique()):
        cluster_bins = clusters[clusters["cluster"] == cluster_id]
        cluster_bed = os.path.join(out_dir, f"cluster_{cluster_id}.bed")
        cluster_bins[["chrom", "start", "end"]].to_csv(cluster_bed, sep="\t", header=False, index=False)
        cluster_total = len(cluster_bins)

        cluster_results = []
        for bed in bed_files:
            if not os.path.exists(bed):
                continue
            annot = os.path.basename(bed)
            p1 = subprocess.Popen(["bedtools", "intersect", "-u", "-a", cluster_bed, "-b", bed], stdout=subprocess.PIPE, text=True)
            p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
            p1.stdout.close()
            overlap_cluster = int((p2.communicate()[0] or "0").strip())

            p1 = subprocess.Popen(["bedtools", "intersect", "-u", "-a", bg_bed, "-b", bed], stdout=subprocess.PIPE, text=True)
            p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
            p1.stdout.close()
            overlap_bg = int((p2.communicate()[0] or "0").strip())
            a = overlap_cluster
            b = cluster_total - overlap_cluster
            c = overlap_bg - overlap_cluster
            d = (bg_total - overlap_bg) - b
            if d < 0:
                d = 0
            odds, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
            cluster_results.append({
                "cluster": cluster_id,
                "annotation": annot,
                "overlap_cluster": a,
                "cluster_total": cluster_total,
                "overlap_background": overlap_bg,
                "background_total": bg_total,
                "odds_ratio": odds,
                "pval": pval,
            })
        if cluster_results:
            pvals = [row["pval"] for row in cluster_results]
            fdrs = bh_fdr(np.array(pvals))
            for row, fdr in zip(cluster_results, fdrs):
                row["fdr"] = fdr
            out_path = os.path.join(out_dir, f"cluster_{cluster_id}.enrichment.tsv")
            pd.DataFrame(cluster_results).to_csv(out_path, sep="\t", index=False)
            all_results.extend(cluster_results)

    if all_results:
        pd.DataFrame(all_results).to_csv(os.path.join(out_dir, "enrichment_summary.tsv"), sep="\t", index=False)
        return "RUN"
    return "NOT_RUN"


def write_summary(path, group, treated, control, n_bins, n_fdr, n_up, n_down, n_clusters, minpts, minsamps, status, reason, enrichment_status):
    header = [
        "group",
        "caller",
        "treated",
        "control",
        "n_bins_tested",
        "n_fdr_pass",
        "n_up",
        "n_down",
        "n_clusters",
        "chosen_minPts",
        "chosen_minSamps",
        "status",
        "reason",
        "enrichment_status",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerow([
            group,
            "NA",
            treated,
            control,
            n_bins,
            n_fdr,
            n_up,
            n_down,
            n_clusters,
            minpts,
            minsamps,
            status,
            reason,
            enrichment_status,
        ])


def write_error_artifact(outdir, message):
    error_path = os.path.join(outdir, "chipbinner.error.txt")
    with open(error_path, "w") as handle:
        handle.write(f"{message}\n")
    return error_path


def main():
    args = parse_args()
    ensure_dir(args.outdir)
    plots_dir = os.path.join(args.outdir, "plots")
    enrichment_dir = os.path.join(args.outdir, "enrichment")

    contrast = [c.strip() for c in args.contrast.split(",") if c.strip()]
    if len(contrast) != 2:
        raise RuntimeError("Contrast must contain treated and control labels")
    treated_label, control_label = contrast

    samples = read_tsv(args.samples)
    if not samples:
        raise RuntimeError("No samples provided")

    samples = [row for row in samples if row.get("condition") in (treated_label, control_label)]
    if not samples:
        raise RuntimeError("No samples match contrast labels")

    sample_ids = [row.get("sample_id", "") for row in samples]

    samplesheet_path = os.path.join(args.outdir, "chipbinner.samplesheet.csv")
    write_samplesheet(samples, samplesheet_path)

    if args.use_input:
        missing_samples = []
        missing_paths = []
        for row in samples:
            bam = row.get("input_bam")
            if not bam or bam in ("NA", "None", "none", ""):
                missing_samples.append(row.get("sample_id", ""))
            elif not os.path.exists(bam):
                missing_paths.append(bam)
        if missing_samples or missing_paths:
            details = []
            if missing_samples:
                details.append(f"samples={','.join(sorted(set(missing_samples)))}")
            if missing_paths:
                details.append(f"paths={','.join(sorted(set(missing_paths)))}")
            reason = "MISSING_INPUT_BAM"
            if details:
                reason = f"{reason}:" + ";".join(details)
            if args.allow_partial:
                error_path = write_error_artifact(args.outdir, reason)
                write_summary(
                    os.path.join(args.outdir, "chipbinner.summary.tsv"),
                    args.group,
                    treated_label,
                    control_label,
                    0,
                    0,
                    0,
                    0,
                    0,
                    "NA",
                    "NA",
                    "SKIP",
                    f"{reason};error={os.path.basename(error_path)}",
                    "NOT_RUN",
                )
                return
            raise RuntimeError(
                f"chipbinner_use_input requested but input BAMs are missing ({'; '.join(details) or 'no input_bam values'})."
            )

    windows_path = os.path.join(args.outdir, "chipbinner.windows.bed")
    try:
        windows_source = choose_windows_file(
            args.windows, args.windows_dir, args.bin_size, args.chrom_sizes, args.blacklist
        )
        if windows_source:
            if args.blacklist and os.path.exists(args.blacklist):
                with tempfile.TemporaryDirectory() as tmpdir:
                    filtered_path = os.path.join(tmpdir, "windows_filtered.bed")
                    run_cmd(["bedtools", "subtract", "-a", windows_source, "-b", args.blacklist], stdout=open(filtered_path, "w"))
                    run_cmd(["sort", "-k1,1", "-k2,2n", filtered_path], stdout=open(windows_path, "w"))
            else:
                run_cmd(["sort", "-k1,1", "-k2,2n", windows_source], stdout=open(windows_path, "w"))
        else:
            generate_windows(args.chrom_sizes, args.bin_size, args.blacklist, windows_path)

        windows_df = load_windows(windows_path)

        counts_raw_path = os.path.join(args.outdir, "chipbinner.bin_counts.tsv")
        counts_df, sample_ids = compute_counts(windows_path, samples, counts_raw_path, use_input=args.use_input)
        write_matrix(counts_df, sample_ids, counts_raw_path)

        norm_matrix, factors_df = normalize_counts(counts_df, sample_ids, samples, args.use_spikein, args.pseudocount)
        norm_path = os.path.join(args.outdir, "chipbinner.normalized_matrix.tsv")
        norm_matrix.insert(0, "bin_id", counts_df["bin_id"])
        norm_matrix.to_csv(norm_path, sep="\t", index=False)

        factors_path = os.path.join(args.outdir, "chipbinner.normalization_factors.tsv")
        factors_df.to_csv(factors_path, sep="\t", index=False)

        ensure_dir(plots_dir)
        plot_pca(norm_matrix.drop(columns=["bin_id"]).values, sample_ids, samples, os.path.join(plots_dir, "chipbinner.pca.png"))
        plot_correlation(norm_matrix.drop(columns=["bin_id"]).values, sample_ids, os.path.join(plots_dir, "chipbinner.correlation.png"))

        n_bins = len(windows_df)
        if n_bins < 2:
            raise RuntimeError("Not enough bins for HDBSCAN clustering (n_bins < 2)")

        def clamp_grid(values, n_bins, min_value):
            vals = sorted({v for v in values if v <= n_bins})
            if not vals:
                fallback = min(n_bins, max(min_value, min(10, n_bins)))
                vals = [fallback]
            return vals

        minpts_list = clamp_grid(parse_grid(args.grid_minpts, [100, 200, 500, 1000]), n_bins, 2)
        minsamps_list = clamp_grid(parse_grid(args.grid_minsamps, [100, 200, 500, 1000]), n_bins, 1)
        grid_dir = os.path.join(args.outdir, "hdbscan_grid")
        grid_summary_path = os.path.join(args.outdir, "chipbinner.hdbscan_grid_summary.tsv")
        chosen, grid_summary_path = hdbscan_grid(norm_matrix.drop(columns=["bin_id"]), windows_df, minpts_list, minsamps_list, grid_dir, grid_summary_path)

        # chosen clusters
        chosen_df = windows_df.copy()
        chosen_df["cluster"] = chosen["labels"]
        clusters_path = os.path.join(args.outdir, "chipbinner.clusters.tsv")
        chosen_df.to_csv(clusters_path, sep="\t", index=False)
        chosen_df.to_csv(os.path.join(args.outdir, "chipbinner.clusters.best.tsv"), sep="\t", index=False)

        # Differential with ROTS
        rots_out = os.path.join(args.outdir, "chipbinner.rots.tsv")
        run_rots(norm_path, samplesheet_path, treated_label, control_label, args.bootstrap, args.k_value, rots_out)
        rots_df = pd.read_csv(rots_out, sep="\t")

        norm_vals = norm_matrix.drop(columns=["bin_id"]).values
        treated_idx = [i for i, row in enumerate(samples) if row.get("condition") == treated_label]
        control_idx = [i for i, row in enumerate(samples) if row.get("condition") == control_label]
        treated_mean = norm_vals[:, treated_idx].mean(axis=1)
        control_mean = norm_vals[:, control_idx].mean(axis=1)
        log2fc = np.log2(treated_mean) - np.log2(control_mean)

        diff_df = pd.DataFrame({
            "chrom": windows_df["chrom"],
            "start": windows_df["start"],
            "end": windows_df["end"],
            "bin_id": windows_df["bin_id"],
            "log2FC": log2fc,
        })
        diff_df = diff_df.merge(rots_df, on="bin_id", how="left")
        diff_df["cluster"] = chosen_df["cluster"].values

        # Standardized cluster views based on cluster mean log2FC (|mean| < 0.25 => stable)
        STABLE_LFC_THRESHOLD = 0.25
        cluster_means = {}
        for cluster_id in sorted(set(chosen_df["cluster"]) - {-1}):
            cluster_means[cluster_id] = float(np.mean(log2fc[chosen_df["cluster"] == cluster_id]))

        def label_two(cluster_id):
            if cluster_id == -1:
                return "noise"
            return "treated_high" if cluster_means.get(cluster_id, 0.0) >= 0 else "control_high"

        def label_three(cluster_id):
            if cluster_id == -1:
                return "noise"
            mean_val = cluster_means.get(cluster_id, 0.0)
            if mean_val >= STABLE_LFC_THRESHOLD:
                return "treated_high"
            if mean_val <= -STABLE_LFC_THRESHOLD:
                return "control_high"
            return "stable"

        cluster_mean_series = chosen_df["cluster"].map(cluster_means).fillna(np.nan)
        clusters_two = chosen_df.copy()
        clusters_two["cluster_label"] = clusters_two["cluster"].map(label_two)
        clusters_two["cluster_mean_log2FC"] = cluster_mean_series
        clusters_two.to_csv(os.path.join(args.outdir, "chipbinner.clusters.2cluster.tsv"), sep="\t", index=False)

        clusters_three = chosen_df.copy()
        clusters_three["cluster_label"] = clusters_three["cluster"].map(label_three)
        clusters_three["cluster_mean_log2FC"] = cluster_mean_series
        clusters_three.to_csv(os.path.join(args.outdir, "chipbinner.clusters.3cluster.tsv"), sep="\t", index=False)

        direction = []
        for _, row in diff_df.iterrows():
            if row.get("FDR") is not None and float(row["FDR"]) <= args.fdr:
                if row["log2FC"] >= args.lfc:
                    direction.append("UP")
                elif row["log2FC"] <= -args.lfc:
                    direction.append("DOWN")
                else:
                    direction.append("NS")
            else:
                direction.append("NS")
        diff_df["direction"] = direction

        diff_path = os.path.join(args.outdir, "chipbinner.differential.tsv")
        diff_df.to_csv(diff_path, sep="\t", index=False)

        n_bins = len(diff_df)
        n_fdr = int((diff_df["FDR"] <= args.fdr).sum()) if "FDR" in diff_df else 0
        n_up = int((diff_df["direction"] == "UP").sum())
        n_down = int((diff_df["direction"] == "DOWN").sum())
        n_clusters = len({lab for lab in chosen["labels"] if lab >= 0})

        enrichment_status = run_enrichment(windows_df, clusters_path, args.functional_db, enrichment_dir)

        summary_path = os.path.join(args.outdir, "chipbinner.summary.tsv")
        write_summary(summary_path, args.group, treated_label, control_label, n_bins, n_fdr, n_up, n_down, n_clusters, chosen["minPts"], chosen["minSamps"], "RUN", "ok", enrichment_status)

    except Exception as exc:
        if args.allow_partial:
            error_path = write_error_artifact(args.outdir, exc)
            write_summary(
                os.path.join(args.outdir, "chipbinner.summary.tsv"),
                args.group,
                treated_label,
                control_label,
                0,
                0,
                0,
                0,
                0,
                "NA",
                "NA",
                "SKIP",
                f"RUNTIME_ERROR:{os.path.basename(error_path)}",
                "NOT_RUN",
            )
        else:
            raise


if __name__ == "__main__":
    main()
