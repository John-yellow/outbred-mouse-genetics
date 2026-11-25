#!/usr/bin/env python3
"""Joint-call the autosomes and X non-PAR in 10-Mb windows using the study settings."""
import argparse
import csv
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

AUTOSOMES = [str(i) for i in range(1, 20)]
X_NONPAR = [(1, 168752754), (169376593, 169476592)]
WINDOW_BP = 10_000_000
ANNOTATE = "FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,INFO/AD"


def read_manifest(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit("Sample manifest is empty")
    ids = [r["sample_id"] for r in rows]
    if len(ids) != len(set(ids)):
        raise SystemExit("Duplicate sample_id in manifest")
    if any(r.get("sex") not in ("M", "F") for r in rows):
        raise SystemExit("Manifest sex values must be M or F")
    return rows


def fai_lengths(path):
    lengths = {}
    for line in path.read_text().splitlines():
        if line:
            col, length, *_ = line.split("\t")
            lengths[col] = int(length)
    return lengths


def make_windows(fai):
    windows = []
    for chrom in AUTOSOMES:
        if chrom not in fai:
            raise SystemExit("Reference .fai is missing chromosome " + chrom)
        for start in range(1, fai[chrom] + 1, WINDOW_BP):
            end = min(fai[chrom], start + WINDOW_BP - 1)
            windows.append((chrom, start, end))
    if "X" not in fai:
        raise SystemExit("Reference .fai is missing chromosome X")
    for low, high in X_NONPAR:
        if high > fai["X"]:
            raise SystemExit("Configured X non-PAR interval exceeds reference length")
        for start in range(low, high + 1, WINDOW_BP):
            end = min(high, start + WINDOW_BP - 1)
            windows.append(("X", start, end))
    return windows


def verify_bcf(bcftools, path, expected_samples):
    """Verify BCF readability, index, non-negative record count, and sample count (bcftools 1.10.2 compatible)."""
    subprocess.run([bcftools, "view", "-h", str(path)], check=True, stdout=subprocess.DEVNULL)
    subprocess.run([bcftools, "index", "-f", str(path)], check=True)
    count_str = subprocess.run([bcftools, "index", "-n", str(path)], check=True,
                               text=True, stdout=subprocess.PIPE).stdout.strip()
    if not count_str.isdigit():
        raise RuntimeError(f"Invalid record count from bcftools index -n for {path}: {count_str!r}")
    samples = subprocess.run([bcftools, "query", "-l", str(path)], check=True,
                             text=True, stdout=subprocess.PIPE).stdout.splitlines()
    if len(samples) != expected_samples:
        raise RuntimeError(f"Sample count mismatch in {path}: {len(samples)} != {expected_samples}")
    return int(count_str)


def run_window(row, bam_list, sex_file, ploidy, ref, outdir, bcftools, call_threads, expected_samples):
    chrom, start, end = row
    label = f"{chrom}_{start}_{end}"
    bcf = outdir / "windows" / f"{label}.bcf"
    logs = outdir / "logs" / label
    logs.mkdir(parents=True, exist_ok=True)
    if bcf.exists() or Path(str(bcf) + ".csi").exists():
        raise RuntimeError(f"Refusing to overwrite existing window output: {bcf}")
    partial = outdir / "windows" / f"{label}.partial.bcf"
    region = f"{chrom}:{start}-{end}"
    mpileup = [bcftools, "mpileup", "-q", "20", "-Q", "20", "-E", "-d", "500",
               "--ff", "UNMAP,SECONDARY,QCFAIL,DUP", "-p", "-m", "3", "-F", "0.25",
               "-a", ANNOTATE, "-f", str(ref), "-r", region, "-b", str(bam_list), "-Ou"]
    call = [bcftools, "call", "-mAv", "-f", "GQ,GP", "--samples-file", str(sex_file),
            "--ploidy-file", str(ploidy), "--threads", str(call_threads), "-Ob", "-o", str(partial)]
    with (logs / "mpileup.stderr.log").open("wb") as e1, (logs / "call.stderr.log").open("wb") as e2:
        p1 = subprocess.Popen(mpileup, stdout=subprocess.PIPE, stderr=e1)
        try:
            p2 = subprocess.Popen(call, stdin=p1.stdout, stdout=subprocess.DEVNULL, stderr=e2)
            p1.stdout.close()
            rc2 = p2.wait()
            rc1 = p1.wait()
        except Exception:
            p1.terminate()
            p1.wait()
            raise
    if rc1 or rc2:
        partial.unlink(missing_ok=True)
        raise RuntimeError(f"Calling failed for {region}: mpileup={rc1}, call={rc2}")
    verify_bcf(bcftools, partial, expected_samples)
    os.replace(partial, bcf)
    os.replace(str(partial) + ".csi", str(bcf) + ".csi")
    return (chrom, start, end, bcf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--bam-dir", type=Path, required=True)
    ap.add_argument("--ploidy", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--call-threads", type=int, default=2)
    args = ap.parse_args()
    for tool in (args.bcftools,):
        subprocess.run([tool, "--version"], check=True, stdout=subprocess.DEVNULL)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run([args.samtools, "faidx", str(args.reference)], check=True)
    manifest = read_manifest(args.manifest)
    bam_paths = [args.bam_dir / f"{r['sample_id']}.marked.bam" for r in manifest]
    for bam in bam_paths:
        bai_options = [Path(str(bam) + ".bai"), bam.with_suffix(".bai")]
        if not bam.is_file() or not any(p.is_file() for p in bai_options):
            raise SystemExit("Missing marked BAM or index: " + str(bam))
    bam_list = args.output_dir / "bams.list"
    bam_list.write_text("".join(str(x.resolve()) + "\n" for x in bam_paths))
    sex_file = args.output_dir / "sex.tsv"
    sex_file.write_text("".join(f"{r['sample_id']}\t{r['sex']}\n" for r in manifest))
    windows = make_windows(fai_lengths(Path(str(args.reference) + ".fai")))
    with (args.output_dir / "windows.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "start", "end", "region"])
        writer.writerows((c, s, e, f"{c}:{s}-{e}") for c, s, e in windows)
    results = {}
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(run_window, row, bam_list, sex_file, args.ploidy,
                               args.reference, args.output_dir, args.bcftools,
                               args.call_threads, len(manifest)): row for row in windows}
        for future in as_completed(futures):
            c, s, e, bcf = future.result()
            results[(c, s, e)] = bcf
            print(f"DONE {c}:{s}-{e}", flush=True)
    ordered = [results[row] for row in windows]
    concat_list = args.output_dir / "window_bcf.list"
    concat_list.write_text("".join(str(p.resolve()) + "\n" for p in ordered))
    final = args.output_dir / "joint_call.raw.bcf"
    subprocess.run([args.bcftools, "concat", "-f", str(concat_list), "-Ob", "-o", str(final)], check=True)
    count = verify_bcf(args.bcftools, final, len(manifest))
    print(f"COMPLETE windows={len(windows)} samples={len(manifest)} records={count} output={final}")


if __name__ == "__main__":
    main()
