#!/usr/bin/env python3
import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path


def detect_cols(header):
    header_lower = [h.lower() for h in header]
    def find(names):
        for name in names:
            if name in header_lower:
                return header[header_lower.index(name)]
        return None
    chr_col = find(["chr", "chrom", "seqnames", "seqname"])
    start_col = find(["start"])
    end_col = find(["end"])
    return chr_col, start_col, end_col


def run(cmd, stdout=None):
    result = subprocess.run(cmd, stdout=stdout, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{result.stderr}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--features", default=None)
    parser.add_argument("--gtf", default=None)
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    with input_path.open() as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader)

    chr_col, start_col, end_col = detect_cols(header)
    if not chr_col or not start_col or not end_col:
        raise SystemExit("Could not detect chr/start/end columns for annotation")

    features_bed = None
    if args.features and Path(args.features).exists():
        features_bed = Path(args.features)
    elif args.gtf and Path(args.gtf).exists():
        raw_bed = Path("features.raw.bed")
        tss_bed = Path("features.tss.bed")
        run(["gtf2bed", args.gtf], stdout=raw_bed.open("w"))
        with raw_bed.open() as in_handle, tss_bed.open("w") as out_handle:
            for line in in_handle:
                if not line.strip():
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) < 6:
                    continue
                chrom, start, end, name, score, strand = fields[:6]
                start = int(start)
                end = int(end)
                if strand == "+":
                    tss_start = start
                else:
                    tss_start = max(end - 1, 0)
                tss_end = tss_start + 1
                out_handle.write("\t".join([chrom, str(tss_start), str(tss_end), name]) + "\n")
        features_bed = tss_bed

    if not features_bed or not features_bed.exists():
        # no annotation possible; copy input to output
        output_path.write_text(input_path.read_text())
        return

    regions_bed = Path("regions.bed")
    with input_path.open() as handle, regions_bed.open("w") as out_handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for idx, row in enumerate(reader):
            out_handle.write("\t".join([
                row[chr_col],
                str(row[start_col]),
                str(row[end_col]),
                str(idx)
            ]) + "\n")

    closest_out = Path("closest.tsv")
    run(["bedtools", "closest", "-d", "-a", str(regions_bed), "-b", str(features_bed)], stdout=closest_out.open("w"))

    annotations = {}
    with closest_out.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 5:
                continue
            row_id = int(fields[3])
            feature_id = fields[7] if len(fields) > 7 else "NA"
            distance = fields[-1]
            annotations[row_id] = (feature_id, distance)

    with input_path.open() as handle, output_path.open("w") as out_handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fieldnames = reader.fieldnames + ["nearest_feature_id", "distance_to_feature"]
        writer = csv.DictWriter(out_handle, delimiter='\t', fieldnames=fieldnames, lineterminator='\n')
        writer.writeheader()
        for idx, row in enumerate(reader):
            feature_id, distance = annotations.get(idx, ("NA", "NA"))
            row["nearest_feature_id"] = feature_id
            row["distance_to_feature"] = distance
            writer.writerow(row)


if __name__ == "__main__":
    main()
