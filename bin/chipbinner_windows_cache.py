#!/usr/bin/env python3
import argparse
import hashlib
import os
import subprocess
import tempfile


def run_cmd(cmd, stdout_path=None):
    if stdout_path:
        with open(stdout_path, "w") as handle:
            result = subprocess.run(cmd, stdout=handle, stderr=subprocess.PIPE, text=True)
    else:
        result = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({' '.join(cmd)}): {result.stderr.strip()}")


def genome_id_from_path(path):
    base = os.path.basename(path)
    if base.endswith(".sizes"):
        base = base[: -len(".sizes")]
    else:
        base = os.path.splitext(base)[0]
    return base or "genome"


def hash_file(path):
    if not path or not os.path.exists(path):
        return "none"
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()[:12]


def resolve_windows_source(windows_dir, genome_id, bin_size, blacklist_hash):
    if not windows_dir:
        return None
    if os.path.isfile(windows_dir):
        return windows_dir
    if not os.path.isdir(windows_dir):
        raise RuntimeError(f"--windows-dir must be a file or directory (got: {windows_dir})")
    preferred = [
        f"windows.{genome_id}.{bin_size}.{blacklist_hash}.bed",
        f"windows.{genome_id}.{bin_size}.bed",
    ]
    for name in preferred:
        candidate = os.path.join(windows_dir, name)
        if os.path.exists(candidate):
            return candidate
    raise RuntimeError(
        "No matching windows BED found in "
        f"{windows_dir}. Expected one of: {', '.join(preferred)}"
    )


def main():
    parser = argparse.ArgumentParser(description="Cache ChIPBinner windows per genome/bin/blacklist")
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--bin-size", type=int, required=True)
    parser.add_argument("--windows-dir")
    parser.add_argument("--blacklist")
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    genome_id = genome_id_from_path(args.chrom_sizes)
    blacklist_hash = hash_file(args.blacklist) if args.blacklist else "none"

    windows_source = resolve_windows_source(args.windows_dir, genome_id, args.bin_size, blacklist_hash)

    out_name = f"windows.{genome_id}.{args.bin_size}.{blacklist_hash}.bed"
    out_path = os.path.join(args.outdir, out_name)

    with tempfile.TemporaryDirectory() as tmpdir:
        base_path = None
        if windows_source:
            base_path = windows_source
        else:
            base_path = os.path.join(tmpdir, "windows_raw.bed")
            run_cmd(["bedtools", "makewindows", "-g", args.chrom_sizes, "-w", str(args.bin_size)], stdout_path=base_path)

        filtered_path = base_path
        if args.blacklist and os.path.exists(args.blacklist):
            filtered_path = os.path.join(tmpdir, "windows_filtered.bed")
            run_cmd(["bedtools", "subtract", "-a", base_path, "-b", args.blacklist], stdout_path=filtered_path)

        run_cmd(["sort", "-k1,1", "-k2,2n", filtered_path], stdout_path=out_path)

    meta_path = os.path.join(args.outdir, "chipbinner_windows_meta.tsv")
    with open(meta_path, "w") as handle:
        handle.write(
            "genome_id\tbin_size\tblacklist_hash\twindows_source\twindows_path\n"
        )
        handle.write(
            f"{genome_id}\t{args.bin_size}\t{blacklist_hash}\t{windows_source or 'generated'}\t{out_name}\n"
        )


if __name__ == "__main__":
    main()
