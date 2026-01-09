#!/usr/bin/env python3
import argparse
import csv
import os


def read_tsv(path):
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for row in reader]


def main():
    parser = argparse.ArgumentParser(description="Build differential summary table")
    parser.add_argument("--design", required=True)
    parser.add_argument("--summaries", nargs="*", default=[])
    parser.add_argument("--manifest-only", action="store_true")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    design_rows = read_tsv(args.design)

    base_rows = {}
    for row in design_rows:
        group = row.get("group") or ""
        caller = row.get("caller") or "NA"
        treated = row.get("treated_condition") or ""
        control = row.get("control_condition") or ""
        status = row.get("status") or ""
        reason = row.get("reason") or ""
        if args.manifest_only and status == "RUN":
            status = "SKIP"
            reason = "manifest_only"

        if caller != "NA":
            key = ("diffbind", group, caller)
            base_rows[key] = {
                "method": "diffbind",
                "group": group,
                "caller": caller,
                "treated": treated,
                "control": control,
                "n_tested": "0",
                "n_fdr_pass": "0",
                "n_up": "0",
                "n_down": "0",
                "n_clusters": "NA",
                "chosen_minPts": "NA",
                "chosen_minSamps": "NA",
                "mode": "NA",
                "status": status,
                "reason": reason,
            }
        else:
            for method in ("chipbinner", "span"):
                method_status = status
                method_reason = reason
                if status == "RUN":
                    if method == "chipbinner" and (row.get("eligible_chipbinner") or "").lower() != "true":
                        method_status = "SKIP"
                        method_reason = "chipbinner_disabled"
                    if method == "span" and (row.get("eligible_span") or "").lower() != "true":
                        method_status = "SKIP"
                        method_reason = "span_disabled"
                key = (method, group, "NA")
                base_rows[key] = {
                    "method": method,
                    "group": group,
                    "caller": "NA",
                    "treated": treated,
                    "control": control,
                    "n_tested": "0",
                    "n_fdr_pass": "0",
                    "n_up": "0",
                    "n_down": "0",
                    "n_clusters": "NA",
                    "chosen_minPts": "NA",
                    "chosen_minSamps": "NA",
                    "mode": "NA",
                    "status": method_status,
                    "reason": method_reason,
                }

    for summary_path in args.summaries:
        if not summary_path or not os.path.exists(summary_path):
            continue
        rows = read_tsv(summary_path)
        if not rows:
            continue
        fname = os.path.basename(summary_path)
        if "diffbind" in fname:
            method = "diffbind"
        elif "chipbinner" in fname:
            method = "chipbinner"
        elif "span" in fname:
            method = "span"
        else:
            continue
        for row in rows:
            group = row.get("group") or ""
            caller = row.get("caller") or (row.get("caller_id") or "NA")
            key = (method, group, caller)
            summary_row = {
                "method": method,
                "group": group,
                "caller": caller,
                "treated": row.get("treated") or row.get("treated_condition") or "",
                "control": row.get("control") or row.get("control_condition") or "",
                "n_tested": row.get("n_tested") or row.get("n_bins_tested") or "0",
                "n_fdr_pass": row.get("n_fdr_pass") or "0",
                "n_up": row.get("n_up") or "0",
                "n_down": row.get("n_down") or "0",
                "n_clusters": row.get("n_clusters") or "NA",
                "chosen_minPts": row.get("chosen_minPts") or "NA",
                "chosen_minSamps": row.get("chosen_minSamps") or "NA",
                "mode": row.get("mode") or "NA",
                "status": row.get("status") or "RUN",
                "reason": row.get("reason") or "ok",
            }
            existing = base_rows.get(key)
            if existing:
                existing_status = (existing.get("status") or "").upper()
                incoming_status = (summary_row.get("status") or "RUN").upper()
                if existing_status not in ("", "RUN") and incoming_status == "RUN":
                    continue
            base_rows[key] = summary_row

    header = [
        "method",
        "group",
        "caller",
        "treated",
        "control",
        "n_tested",
        "n_fdr_pass",
        "n_up",
        "n_down",
        "n_clusters",
        "chosen_minPts",
        "chosen_minSamps",
        "mode",
        "status",
        "reason",
    ]

    rows_sorted = sorted(base_rows.values(), key=lambda r: (r["method"], r["group"], r["caller"]))

    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_sorted)


if __name__ == "__main__":
    main()
