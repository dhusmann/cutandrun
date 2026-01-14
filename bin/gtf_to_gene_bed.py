#!/usr/bin/env python3
import argparse
import gzip


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_attributes(attr_str):
    attrs = {}
    for raw in attr_str.strip().strip(";").split(";"):
        item = raw.strip()
        if not item:
            continue
        if "=" in item and " " not in item:
            key, val = item.split("=", 1)
        else:
            parts = item.split(" ", 1)
            if len(parts) < 2:
                continue
            key, val = parts[0], parts[1]
        val = val.strip().strip('"')
        attrs[key] = val
    return attrs


def pick_attr(attrs, keys):
    for key in keys:
        val = attrs.get(key)
        if val:
            return val
    return "NA"


def main():
    parser = argparse.ArgumentParser(description="Create gene BED from GTF")
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--feature", default="gene")
    args = parser.parse_args()

    with open_text(args.gtf) as handle, open(args.out, "w") as out:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _source, feature, start, end, _score, _strand, _frame, attrs_str = fields[:9]
            if feature != args.feature:
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            if end_i < start_i:
                start_i, end_i = end_i, start_i
            start0 = max(start_i - 1, 0)

            attrs = parse_attributes(attrs_str)
            gene_id = pick_attr(attrs, ["gene_id", "ID", "gene", "transcript_id"])
            gene_name = pick_attr(attrs, ["gene_name", "Name", "gene"])
            if gene_id == "NA" and gene_name != "NA":
                gene_id = gene_name

            out.write(f"{chrom}\t{start0}\t{end_i}\t{gene_id}\t{gene_name}\n")


if __name__ == "__main__":
    main()
