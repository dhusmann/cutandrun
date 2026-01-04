#!/usr/bin/env python3
import argparse
import json


def main():
    parser = argparse.ArgumentParser(description="Write records JSON to TSV")
    parser.add_argument("--records", required=True, help="JSON string of records")
    parser.add_argument("--header", required=True, help="Tab-separated header columns")
    parser.add_argument("--out", required=True, help="Output TSV path")
    args = parser.parse_args()

    records = json.loads(args.records)
    header = args.header.split("\t")

    with open(args.out, "w") as handle:
        handle.write("\t".join(header) + "\n")
        for record in records:
            row = [str(record.get(col, "")) for col in header]
            handle.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
