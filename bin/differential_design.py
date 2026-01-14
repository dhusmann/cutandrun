#!/usr/bin/env python3
import argparse
import csv
import sys


def read_tsv(path):
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for row in reader]


def parse_list(value):
    if not value:
        return None
    items = [x.strip() for x in value.split(",") if x.strip()]
    return set(items) if items else None


def is_span_caller(name):
    if not name:
        return False
    lowered = str(name).lower()
    return lowered.startswith("span") or lowered.startswith("omnipeak")


def main():
    parser = argparse.ArgumentParser(description="Build differential design manifest")
    parser.add_argument("--samples", required=True)
    parser.add_argument("--peaks", required=True)
    parser.add_argument("--contrast", required=True)
    parser.add_argument("--min-replicates", type=int, required=True)
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument("--run-diffbind", action="store_true")
    parser.add_argument("--run-chipbinner", action="store_true")
    parser.add_argument("--run-span", action="store_true")
    parser.add_argument("--groups")
    parser.add_argument("--callers")
    parser.add_argument("--out", default="differential_manifest.design.tsv")
    args = parser.parse_args()

    contrast = [x.strip() for x in args.contrast.split(",") if x.strip()]
    if len(contrast) != 2:
        print("ERROR: --differential_contrast must contain exactly two comma-separated labels.", file=sys.stderr)
        sys.exit(1)
    treated, control = contrast

    samples = read_tsv(args.samples)
    peaks = read_tsv(args.peaks)

    if not samples:
        print("ERROR: differential manifest has no samples.", file=sys.stderr)
        sys.exit(1)

    conditions_present = {row.get("condition") for row in samples if row.get("condition")}
    if treated not in conditions_present or control not in conditions_present:
        print(
            "ERROR: differential_contrast labels not found in samples: "
            f"treated={treated} control={control} present={sorted(conditions_present)}",
            file=sys.stderr,
        )
        sys.exit(1)

    group_allow = parse_list(args.groups)
    caller_allow = parse_list(args.callers)

    # Index samples by group
    samples_by_group = {}
    for row in samples:
        group = row.get("group") or ""
        samples_by_group.setdefault(group, []).append(row)

    # Index peaks by (group, caller, sample_id)
    peaks_by_group_caller = {}
    for row in peaks:
        group = row.get("group") or ""
        caller = row.get("caller") or ""
        sample_id = row.get("sample_id") or ""
        key = (group, caller)
        peaks_by_group_caller.setdefault(key, set()).add(sample_id)

    all_callers = sorted({row.get("caller") or "" for row in peaks})

    span_callers_present = True
    if args.run_span:
        span_callers = [caller for caller in all_callers if is_span_caller(caller)]
        span_callers_present = bool(span_callers)
        if not span_callers_present:
            print(
                "WARNING: run_span_diff requested but no SPAN/OmniPeak peaks found in manifest; "
                "SPAN differential will be skipped.",
                file=sys.stderr,
            )

    design_rows = []
    fail_required = False

    for group in sorted(samples_by_group.keys()):
        group_samples = samples_by_group[group]
        conditions = sorted({row.get("condition") for row in group_samples if row.get("condition")})
        ignored = [c for c in conditions if c not in contrast]
        ignored_conditions = ",".join(ignored) if ignored else "NA"

        treated_samples = [row for row in group_samples if row.get("condition") == treated]
        control_samples = [row for row in group_samples if row.get("condition") == control]

        treated_reps = {row.get("replicate") for row in treated_samples if row.get("replicate")}
        control_reps = {row.get("replicate") for row in control_samples if row.get("replicate")}

        n_treated = len(treated_reps)
        n_control = len(control_reps)

        group_filtered = group_allow is not None and group not in group_allow

        base_status = "RUN"
        base_reason = "ok"
        if group_filtered:
            base_status = "SKIP"
            base_reason = "group_filtered"
        elif not treated_samples or not control_samples:
            base_status = "SKIP"
            base_reason = "missing_condition"
        elif n_treated < args.min_replicates or n_control < args.min_replicates:
            if args.allow_partial:
                base_status = "SKIP"
                base_reason = "insufficient_replicates"
            else:
                base_status = "FAIL"
                base_reason = "insufficient_replicates"

        group_eligible = base_status == "RUN"

        # Group-level row (chipbinner/span)
        if args.run_chipbinner or args.run_span:
            status = base_status
            reason = base_reason
            eligible_chipbinner = args.run_chipbinner and group_eligible
            eligible_span = args.run_span and span_callers_present and group_eligible
            if not (args.run_chipbinner or args.run_span):
                status = "SKIP"
                reason = "methods_disabled"
        else:
            status = "SKIP"
            reason = "methods_disabled"
            eligible_chipbinner = False
            eligible_span = False

        design_rows.append({
            "group": group,
            "caller": "NA",
            "treated_condition": treated,
            "control_condition": control,
            "n_treated": str(n_treated),
            "n_control": str(n_control),
            "eligible_diffbind": "false",
            "eligible_chipbinner": "true" if eligible_chipbinner else "false",
            "eligible_span": "true" if eligible_span else "false",
            "status": status,
            "reason": reason,
            "ignored_conditions": ignored_conditions,
        })

        if status == "FAIL" and (args.run_chipbinner or args.run_span) and not args.allow_partial:
            fail_required = True

        # Diffbind rows per caller
        for caller in all_callers:
            if caller_allow is not None and caller not in caller_allow:
                status = "SKIP"
                reason = "caller_filtered"
                eligible_diffbind = False
            elif not args.run_diffbind:
                status = "SKIP"
                reason = "diffbind_disabled"
                eligible_diffbind = False
            else:
                status = base_status
                reason = base_reason
                eligible_diffbind = group_eligible

                if status == "RUN":
                    key = (group, caller)
                    caller_samples = peaks_by_group_caller.get(key, set())
                    needed_ids = {row.get("sample_id") for row in treated_samples + control_samples}
                    missing = sorted([sid for sid in needed_ids if sid not in caller_samples])
                    if missing:
                        if args.allow_partial:
                            status = "SKIP"
                            reason = "missing_peaks"
                        else:
                            status = "FAIL"
                            reason = "missing_peaks"
                        eligible_diffbind = False

            design_rows.append({
                "group": group,
                "caller": caller,
                "treated_condition": treated,
                "control_condition": control,
                "n_treated": str(n_treated),
                "n_control": str(n_control),
                "eligible_diffbind": "true" if eligible_diffbind else "false",
                "eligible_chipbinner": "true" if (args.run_chipbinner and group_eligible) else "false",
                "eligible_span": "true"
                if (args.run_span and span_callers_present and group_eligible)
                else "false",
                "status": status,
                "reason": reason,
                "ignored_conditions": ignored_conditions,
            })

            if status == "FAIL" and args.run_diffbind and not args.allow_partial:
                fail_required = True

    # Sort rows deterministically
    design_rows.sort(key=lambda r: (r["group"], r["caller"], r["treated_condition"], r["control_condition"]))

    header = [
        "group",
        "caller",
        "treated_condition",
        "control_condition",
        "n_treated",
        "n_control",
        "eligible_diffbind",
        "eligible_chipbinner",
        "eligible_span",
        "status",
        "reason",
        "ignored_conditions",
    ]

    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(design_rows)

    if fail_required:
        print("ERROR: differential design has FAIL entries and --differential_allow_partial is false.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
