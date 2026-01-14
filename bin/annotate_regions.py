#!/usr/bin/env python3
import argparse
import csv
import os
import subprocess
import tempfile


def read_tsv(path):
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = [row for row in reader]
    return rows


def is_header(row):
    if len(row) < 3:
        return True
    try:
        int(row[1])
        return False
    except Exception:
        return True


def main():
    parser = argparse.ArgumentParser(description="Annotate regions with nearest gene")
    parser.add_argument("--regions", required=True)
    parser.add_argument("--gene-bed", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    rows = read_tsv(args.regions)
    if not rows:
        with open(args.out, "w") as handle:
            handle.write("")
        return

    header = None
    data_rows = rows
    if is_header(rows[0]):
        header = rows[0]
        data_rows = rows[1:]

    annotation_cols = ["nearest_feature_id", "nearest_gene_name", "distance_to_feature"]

    if not data_rows:
        with open(args.out, "w") as handle:
            if header:
                handle.write("\t".join(header + annotation_cols) + "\n")
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        bed_path = os.path.join(tmpdir, "regions.bed")
        with open(bed_path, "w") as handle:
            for idx, row in enumerate(data_rows):
                handle.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{idx}\n")

        cmd = [
            "bedtools",
            "closest",
            "-d",
            "-a",
            bed_path,
            "-b",
            args.gene_bed,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        annotations = {}
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            parts = line.split("\t")
            idx = int(parts[3])
            a_cols = 4
            distance = parts[-1] if parts else "NA"
            b_cols = len(parts) - a_cols - 1
            gene_id = "NA"
            gene_name = "NA"
            if b_cols >= 4:
                gene_id = parts[a_cols + 3]
            if b_cols >= 5:
                gene_name = parts[a_cols + 4]
            if gene_id in [".", ""]:
                gene_id = "NA"
            if gene_name in [".", ""]:
                gene_name = "NA"
            annotations[idx] = (gene_id, gene_name, distance)

    with open(args.out, "w") as handle:
        if header:
            handle.write("\t".join(header + annotation_cols) + "\n")
        for idx, row in enumerate(data_rows):
            gene_id, gene_name, distance = annotations.get(idx, ("NA", "NA", "NA"))
            handle.write("\t".join(row + [gene_id, gene_name, distance]) + "\n")


if __name__ == "__main__":
    main()
