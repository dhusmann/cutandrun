#!/usr/bin/env python3
import argparse
import itertools
import sys

import numpy as np
import pandas as pd


def parse_list(value):
    return [int(x.strip()) for x in str(value).split(',') if x.strip()]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--min_cluster_size", required=True)
    parser.add_argument("--min_samples", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--clusters", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.matrix, sep='\t')
    if df.shape[1] <= 3:
        raise SystemExit("Matrix missing sample columns")

    coords = df.iloc[:, :3]
    matrix = df.iloc[:, 3:].values

    try:
        import hdbscan
    except Exception:
        # fallback: single cluster
        summary = pd.DataFrame([
            {
                "min_cluster_size": 0,
                "min_samples": 0,
                "n_clusters": 1,
                "fraction_assigned": 1.0,
                "mean_persistence": 1.0,
                "selected": True,
            }
        ])
        summary.to_csv(args.summary, sep='\t', index=False)
        clusters = coords.copy()
        clusters["cluster"] = 0
        clusters.to_csv(args.clusters, sep='\t', index=False)
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

    summary = pd.DataFrame(results)
    if best is not None:
        summary.loc[(summary.min_cluster_size == best[0]) & (summary.min_samples == best[1]), "selected"] = True
    summary.to_csv(args.summary, sep='\t', index=False)

    clusters = coords.copy()
    clusters["cluster"] = best_labels if best_labels is not None else 0
    clusters.to_csv(args.clusters, sep='\t', index=False)


if __name__ == "__main__":
    main()
