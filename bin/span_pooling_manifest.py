#!/usr/bin/env python3
import argparse
import csv
import os


def read_rows(path):
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for row in reader]


def main():
    parser = argparse.ArgumentParser(description="Merge SPAN pooling manifests")
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = []
    for path in args.inputs:
        if not path or not os.path.exists(path):
            continue
        rows.extend(read_rows(path))

    if not rows:
        return

    rows.sort(key=lambda r: (r.get("group") or "", r.get("condition") or "", r.get("pooled_bam") or ""))

    header = ["group", "condition", "pooled_bam", "source_bams"]
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
