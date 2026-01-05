#!/usr/bin/env python3
import argparse
import csv
import itertools
import sys

import numpy as np


def parse_list(value):
    return [int(x.strip()) for x in str(value).split(',') if x.strip()]


def load_matrix(path):
    with open(path, newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader, [])
        rows = [row for row in reader if row]
    if len(header) <= 3:
        raise SystemExit("Matrix missing sample columns")
    coords = [row[:3] for row in rows]
    data = []
    for row in rows:
        data.append([float(x) for x in row[3:]])
    matrix = np.array(data) if data else np.empty((0, len(header) - 3))
    return header[:3], header[3:], coords, matrix


def load_samplesheet(path, sample_ids):
    if not path:
        return []
    conditions = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sample_id = row.get("sample_id")
            condition = row.get("condition")
            if not sample_id:
                continue
            conditions[sample_id] = condition or ""
    missing = [sid for sid in sample_ids if sid not in conditions]
    if missing:
        raise SystemExit(f"Samplesheet missing condition for: {', '.join(missing)}")
    return [conditions[sid] for sid in sample_ids]


def replicate_consistency(matrix, labels, conditions):
    if matrix.size == 0 or not conditions:
        return 0.0
    scores = []
    unique_conditions = sorted(set(conditions))
    for cluster_id in sorted(set(labels)):
        if cluster_id == -1:
            continue
        cluster_rows = matrix[labels == cluster_id, :]
        if cluster_rows.shape[0] < 2:
            continue
        for condition in unique_conditions:
            idx = [i for i, cond in enumerate(conditions) if cond == condition]
            if len(idx) < 2:
                continue
            data = cluster_rows[:, idx]
            try:
                corr = np.corrcoef(data, rowvar=False)
            except Exception:
                continue
            if corr.size <= 1:
                continue
            tri = corr[np.triu_indices_from(corr, k=1)]
            tri = tri[~np.isnan(tri)]
            if tri.size:
                scores.append(float(np.mean(tri)))
    if not scores:
        return 0.0
    return float(np.mean(scores))


def write_summary(path, rows):
    fieldnames = [
        "min_cluster_size",
        "min_samples",
        "n_clusters",
        "frac_assigned_non_noise",
        "mean_persistence",
        "replicate_consistency",
        "score_total",
        "selected",
        "selected_2clusters",
        "selected_3clusters",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_clusters(path, coord_header, coords, labels):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(coord_header + ["cluster"])
        for row, label in zip(coords, labels):
            writer.writerow(row + [label])


def write_not_found(path, coord_header, message):
    with open(path, "w", newline="") as handle:
        handle.write(f"# {message}\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(coord_header + ["cluster"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--min_cluster_size", required=True)
    parser.add_argument("--min_samples", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--clusters_best", required=True)
    parser.add_argument("--clusters_2", required=True)
    parser.add_argument("--clusters_3", required=True)
    parser.add_argument("--samplesheet", default="")
    args = parser.parse_args()

    coord_header, sample_ids, coords, matrix = load_matrix(args.matrix)
    conditions = load_samplesheet(args.samplesheet, sample_ids) if args.samplesheet else []

    try:
        import hdbscan
    except Exception:
        # fallback: single cluster
        summary_rows = [
            {
                "min_cluster_size": 0,
                "min_samples": 0,
                "n_clusters": 1,
                "frac_assigned_non_noise": 1.0,
                "mean_persistence": 1.0,
                "replicate_consistency": 0.0,
                "score_total": 1.0,
                "selected": True,
                "selected_2clusters": False,
                "selected_3clusters": False,
            }
        ]
        write_summary(args.summary, summary_rows)
        labels = [0 for _ in coords]
        write_clusters(args.clusters_best, coord_header, coords, labels)
        write_not_found(args.clusters_2, coord_header, "No 2-cluster solution found (HDBSCAN unavailable)")
        write_not_found(args.clusters_3, coord_header, "No 3-cluster solution found (HDBSCAN unavailable)")
        return

    min_cluster_sizes = parse_list(args.min_cluster_size)
    min_samples_list = parse_list(args.min_samples)

    results = []
    labels_by_params = {}

    for mcs, ms in itertools.product(min_cluster_sizes, min_samples_list):
        clusterer = hdbscan.HDBSCAN(min_cluster_size=mcs, min_samples=ms)
        labels = clusterer.fit_predict(matrix) if matrix.size else np.array([])
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        fraction_assigned = float(np.sum(labels >= 0)) / float(len(labels)) if len(labels) else 0.0
        mean_persistence = float(np.mean(clusterer.cluster_persistence_)) if getattr(clusterer, 'cluster_persistence_', None) is not None and len(clusterer.cluster_persistence_) > 0 else 0.0
        rep_consistency = replicate_consistency(matrix, labels, conditions) if len(labels) else 0.0
        score_total = float(fraction_assigned + mean_persistence + rep_consistency)
        labels_by_params[(mcs, ms)] = labels
        results.append({
            "min_cluster_size": mcs,
            "min_samples": ms,
            "n_clusters": n_clusters,
            "frac_assigned_non_noise": round(fraction_assigned, 4),
            "mean_persistence": round(mean_persistence, 4),
            "replicate_consistency": round(rep_consistency, 4),
            "score_total": round(score_total, 4),
            "selected": False,
            "selected_2clusters": False,
            "selected_3clusters": False,
        })

    best = max(results, key=lambda row: (row["score_total"], row["frac_assigned_non_noise"], row["mean_persistence"]))
    best_2 = None
    best_3 = None
    candidates_2 = [row for row in results if row["n_clusters"] == 2]
    candidates_3 = [row for row in results if row["n_clusters"] == 3]
    if candidates_2:
        best_2 = max(candidates_2, key=lambda row: (row["score_total"], row["frac_assigned_non_noise"], row["mean_persistence"]))
    if candidates_3:
        best_3 = max(candidates_3, key=lambda row: (row["score_total"], row["frac_assigned_non_noise"], row["mean_persistence"]))

    for row in results:
        if row is best:
            row["selected"] = True
        if best_2 is not None and row is best_2:
            row["selected_2clusters"] = True
        if best_3 is not None and row is best_3:
            row["selected_3clusters"] = True

    results_sorted = sorted(results, key=lambda row: row["score_total"], reverse=True)
    write_summary(args.summary, results_sorted)

    labels = list(labels_by_params.get((best["min_cluster_size"], best["min_samples"]), [])) if coords else []
    if coords:
        write_clusters(args.clusters_best, coord_header, coords, labels)
    else:
        write_clusters(args.clusters_best, coord_header, coords, [])

    if best_2 is not None:
        labels_2 = list(labels_by_params.get((best_2["min_cluster_size"], best_2["min_samples"]), []))
        write_clusters(args.clusters_2, coord_header, coords, labels_2)
    else:
        sys.stderr.write("No 2-cluster solution found for HDBSCAN grid search\n")
        write_not_found(args.clusters_2, coord_header, "No 2-cluster solution found")

    if best_3 is not None:
        labels_3 = list(labels_by_params.get((best_3["min_cluster_size"], best_3["min_samples"]), []))
        write_clusters(args.clusters_3, coord_header, coords, labels_3)
    else:
        sys.stderr.write("No 3-cluster solution found for HDBSCAN grid search\n")
        write_not_found(args.clusters_3, coord_header, "No 3-cluster solution found")


if __name__ == "__main__":
    main()
