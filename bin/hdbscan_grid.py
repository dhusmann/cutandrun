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
    return header[:3], coords, matrix


def write_summary(path, rows):
    fieldnames = [
        "min_cluster_size",
        "min_samples",
        "n_clusters",
        "fraction_assigned",
        "mean_persistence",
        "selected",
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--min_cluster_size", required=True)
    parser.add_argument("--min_samples", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--clusters", required=True)
    args = parser.parse_args()

    coord_header, coords, matrix = load_matrix(args.matrix)

    try:
        import hdbscan
    except Exception:
        # fallback: single cluster
        summary_rows = [
            {
                "min_cluster_size": 0,
                "min_samples": 0,
                "n_clusters": 1,
                "fraction_assigned": 1.0,
                "mean_persistence": 1.0,
                "selected": True,
            }
        ]
        write_summary(args.summary, summary_rows)
        labels = [0 for _ in coords]
        write_clusters(args.clusters, coord_header, coords, labels)
        return

    min_cluster_sizes = parse_list(args.min_cluster_size)
    min_samples_list = parse_list(args.min_samples)

    results = []
    best = None
    best_score = None
    best_labels = None

    for mcs, ms in itertools.product(min_cluster_sizes, min_samples_list):
        clusterer = hdbscan.HDBSCAN(min_cluster_size=mcs, min_samples=ms)
        labels = clusterer.fit_predict(matrix)
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        fraction_assigned = float(np.sum(labels >= 0)) / float(len(labels)) if len(labels) else 0.0
        mean_persistence = float(np.mean(clusterer.cluster_persistence_)) if getattr(clusterer, 'cluster_persistence_', None) is not None and len(clusterer.cluster_persistence_) > 0 else 0.0
        score = (fraction_assigned, n_clusters, mean_persistence)
        if best_score is None or score > best_score:
            best_score = score
            best = (mcs, ms)
            best_labels = labels
        results.append({
            "min_cluster_size": mcs,
            "min_samples": ms,
            "n_clusters": n_clusters,
            "fraction_assigned": round(fraction_assigned, 4),
            "mean_persistence": round(mean_persistence, 4),
            "selected": False,
        })

    for row in results:
        if best is not None and row["min_cluster_size"] == best[0] and row["min_samples"] == best[1]:
            row["selected"] = True
    write_summary(args.summary, results)

    labels = list(best_labels) if best_labels is not None else [0 for _ in coords]
    write_clusters(args.clusters, coord_header, coords, labels)


if __name__ == "__main__":
    main()
