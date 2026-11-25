#!/usr/bin/env python3
"""Validate BCF-to-VCF.gz conversion integrity, sample order, genotypes, X ploidy, and CSI index access."""
import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

PAR_START = 168752755
PAR_END = 169376592


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def query_count(bcftools, path):
    p = subprocess.run(
        [bcftools, "index", "-n", str(path)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return int(p.stdout.strip())


def samples(bcftools, path):
    p = subprocess.run(
        [bcftools, "query", "-l", str(path)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return p.stdout.splitlines()


def load_sex_map(path, expected_samples):
    sex_map = {}
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        sid, sex = parts[0], parts[1]
        if sex in ("M", "F"):
            sex_map[sid] = sex
    missing = [s for s in expected_samples if s not in sex_map]
    if missing:
        raise RuntimeError(f"Sex map is missing {len(missing)} samples, e.g. {missing[:3]}")
    return [sex_map[s] for s in expected_samples]


def check_x_ploidy(line_bytes, sample_sexes):
    cols = line_bytes.rstrip(b"\r\n").decode("utf-8").split("\t")
    pos = int(cols[1])
    if PAR_START <= pos <= PAR_END:
        raise RuntimeError(f"PAR site present in X non-PAR callset at X:{pos}")
    fmt_parts = cols[8].split(":")
    try:
        gt_idx = fmt_parts.index("GT")
    except ValueError as exc:
        raise RuntimeError(f"Missing GT in FORMAT at X:{pos}: {cols[8]}") from exc
    sample_cols = cols[9:]
    if len(sample_cols) != len(sample_sexes):
        raise RuntimeError(f"Sample column count mismatch at X:{pos}: {len(sample_cols)} != {len(sample_sexes)}")
    for idx, (scol, sex) in enumerate(zip(sample_cols, sample_sexes)):
        gt = scol.split(":")[gt_idx]
        if sex == "M":
            if gt not in ("0", "1", "."):
                raise RuntimeError(f"Invalid male X genotype at X:{pos} sample_idx={idx}: {gt!r}")
        else:
            alleles = gt.replace("|", "/").split("/")
            if len(alleles) != 2 or any(a not in ("0", "1", ".") for a in alleles):
                raise RuntimeError(f"Invalid female X genotype at X:{pos} sample_idx={idx}: {gt!r}")


def compare_streams(bcftools, source, target, sample_sexes=None):
    a = subprocess.Popen(
        [bcftools, "view", "--no-version", "-H", str(source)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    b = subprocess.Popen(
        [bcftools, "view", "--no-version", "-H", str(target)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    count = 0
    x_checked = 0
    first = None
    while True:
        x = a.stdout.readline()
        y = b.stdout.readline()
        if not x and not y:
            break
        count += 1
        if x != y:
            first = (
                count,
                x[:200].decode(errors="replace").rstrip(),
                y[:200].decode(errors="replace").rstrip(),
            )
            a.terminate()
            b.terminate()
            break
        if sample_sexes is not None and x.startswith(b"X\t"):
            check_x_ploidy(x, sample_sexes)
            x_checked += 1
    a.stdout.close()
    b.stdout.close()
    ae = a.stderr.read().decode(errors="replace")
    be = b.stderr.read().decode(errors="replace")
    ar = a.wait()
    br = b.wait()
    if first:
        raise RuntimeError(
            f"BCF/VCF full record/genotype mismatch at row {first[0]}: {first[1]} <> {first[2]}"
        )
    if ar or br:
        raise RuntimeError(f"bcftools view -H failed: source={ar} {ae}; target={br} {be}")
    return count, x_checked


def first_pos_on_chrom(bcftools, path, chrom):
    p = subprocess.Popen(
        [bcftools, "query", "-f", "%POS\n", "-r", chrom, str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    first_line = p.stdout.readline().strip()
    p.terminate()
    p.stdout.close()
    p.stderr.close()
    p.wait()
    if not first_line:
        raise RuntimeError(f"No records found on chromosome {chrom} in {path}")
    return int(first_line)


def check_csi_region(bcftools, src, vcf, chrom):
    pos = first_pos_on_chrom(bcftools, src, chrom)
    region = f"{chrom}:{pos}-{pos + 1000}"
    src_out = subprocess.run(
        [bcftools, "view", "--no-version", "-H", "-r", region, str(src)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    vcf_out = subprocess.run(
        [bcftools, "view", "--no-version", "-H", "-r", region, str(vcf)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    src_lines = src_out.splitlines()
    vcf_lines = vcf_out.splitlines()
    if len(vcf_lines) < 1:
        raise RuntimeError(f"CSI region query {region} returned 0 records in {vcf}")
    if src_lines != vcf_lines:
        raise RuntimeError(f"CSI region query {region} mismatch between {src} and {vcf}")
    return {
        "region": region,
        "source_records": len(src_lines),
        "vcf_records": len(vcf_lines),
        "exact_match": True,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-bcf", type=Path, required=True)
    ap.add_argument("--vcf-gz", type=Path, required=True)
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--sex-map", type=Path, default=None, help="Optional TSV (sample_id<TAB>M|F) for X ploidy validation")
    ap.add_argument("--report", type=Path, required=True)
    args = ap.parse_args()

    src = args.source_bcf.resolve()
    vcf = args.vcf_gz.resolve()
    if not src.is_file() or not vcf.is_file():
        raise SystemExit("Source BCF or VCF.gz is missing")
    csi = Path(str(vcf) + ".csi")
    if not csi.is_file():
        raise SystemExit("CSI index is missing: " + str(csi))

    source_samples = samples(args.bcftools, src)
    vcf_samples = samples(args.bcftools, vcf)
    if source_samples != vcf_samples:
        raise RuntimeError("Sample order differs between source BCF and VCF")

    sample_sexes = None
    if args.sex_map is not None:
        sample_sexes = load_sex_map(args.sex_map.resolve(), source_samples)

    source_n = query_count(args.bcftools, src)
    vcf_n = query_count(args.bcftools, vcf)
    if source_n != vcf_n:
        raise RuntimeError(f"Indexed site counts differ: {source_n} != {vcf_n}")

    headers = subprocess.run(
        [args.bcftools, "view", "-h", str(vcf)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    ).stdout
    if "##INFO=<ID=CSQ," not in headers:
        raise RuntimeError("VCF header is missing INFO/CSQ declaration")

    compared, x_checked = compare_streams(args.bcftools, src, vcf, sample_sexes=sample_sexes)
    if compared != source_n:
        raise RuntimeError(f"Streamed count differs from index count: {compared} != {source_n}")

    csi_checks = {
        "1": check_csi_region(args.bcftools, src, vcf, "1"),
        "X": check_csi_region(args.bcftools, src, vcf, "X"),
    }

    report = {
        "status": "PASS",
        "validated_at": datetime.now(timezone.utc).isoformat(),
        "source_bcf": str(src),
        "vcf_gz": str(vcf),
        "csi": str(csi),
        "source_bytes": src.stat().st_size,
        "vcf_gz_bytes": vcf.stat().st_size,
        "csi_bytes": csi.stat().st_size,
        "sample_count": len(source_samples),
        "sample_order": "identical",
        "source_indexed_record_count": source_n,
        "vcf_indexed_record_count": vcf_n,
        "streamed_rows_compared": compared,
        "full_record_and_genotype_mismatches": 0,
        "compared_fields": [
            "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "ALL_SAMPLES"
        ],
        "csq_header": "present",
        "csi_region_checks": csi_checks,
        "sha256_vcf_gz": sha256(vcf),
        "sha256_csi": sha256(csi),
    }
    if sample_sexes is not None:
        report["x_ploidy_checks"] = {
            "sex_map": str(args.sex_map.resolve()),
            "records_checked": x_checked,
            "par_records": 0,
            "ploidy_errors": 0,
        }

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
