#!/usr/bin/env python3
"""Streaming, sample-level DP/GQ mask with ploidy-preserving missing GTs."""
import argparse
import collections
import sys

PAR_X = ((168752755, 169376592),)
F_MISSING_HEADER = '##INFO=<ID=F_MISSING,Number=1,Type=Float,Description="Fraction of samples with missing GT after per-sample DP/GQ masking">'


def load_sexes(path):
    sexes = []
    seen = set()
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 2 or fields[1] not in ('M', 'F'):
                raise ValueError('sex map must contain sample<TAB>M/F: ' + line.rstrip())
            if fields[0] in seen:
                raise ValueError('duplicate sample in sex map: ' + fields[0])
            seen.add(fields[0])
            sexes.append((fields[0], fields[1]))
    return sexes


def expected_ploidy(chrom, pos, sex):
    if chrom == 'X':
        if any(start <= pos <= end for start, end in PAR_X):
            raise ValueError('PAR record entered non-PAR callset at X:' + str(pos))
        return 1 if sex == 'M' else 2
    if chrom in set(str(i) for i in range(1, 20)):
        return 2
    raise ValueError('unexpected contig in filtering input: ' + chrom)


def canonical_missing(ploidy):
    return '.' if ploidy == 1 else './.'


def parse_int(value, label, chrom, pos, sample):
    if value in ('', '.'):
        return None
    try:
        return int(value)
    except ValueError:
        raise ValueError('{} is non-integer at {}:{} sample {}: {}'.format(label, chrom, pos, sample, value))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sex-map', required=True)
    args = parser.parse_args()
    sexes = load_sexes(args.sex_map)
    sample_names = None
    total = low_dp = low_gq = missing_metric = originally_missing = masked = 0
    chrom_counts = collections.Counter()
    fmissing_sum = 0.0
    wrote_fmissing_header = False

    for line in sys.stdin:
        if line.startswith('##INFO=<ID=F_MISSING,'):
            # Caller INFO is not trusted after genotype edits.
            continue
        if line.startswith('#CHROM\t'):
            fields = line.rstrip('\n').split('\t')
            sample_names = fields[9:]
            expected_names = [name for name, _ in sexes]
            if sample_names != expected_names:
                raise SystemExit('VCF samples/order differ from frozen sex map')
            sys.stdout.write(F_MISSING_HEADER + '\n')
            wrote_fmissing_header = True
            sys.stdout.write(line)
            continue
        if line.startswith('#'):
            sys.stdout.write(line)
            continue
        if sample_names is None or not wrote_fmissing_header:
            raise SystemExit('VCF header missing #CHROM before records')
        fields = line.rstrip('\n').split('\t')
        if len(fields) != 9 + len(sexes):
            raise SystemExit('VCF record has unexpected column count at ' + ':'.join(fields[:2]))
        chrom, pos = fields[0], int(fields[1])
        fmt = fields[8].split(':')
        if not all(key in fmt for key in ('GT', 'DP', 'GQ')):
            raise SystemExit('GT, DP, or GQ missing from FORMAT at {}:{}'.format(chrom, pos))
        gt_idx, dp_idx, gq_idx = (fmt.index(x) for x in ('GT', 'DP', 'GQ'))
        info_parts = [] if fields[7] in ('', '.') else [x for x in fields[7].split(';') if not x.startswith('F_MISSING=')]
        miss_site = 0
        for sample_i, (sample, sex) in enumerate(sexes):
            ploidy = expected_ploidy(chrom, pos, sex)
            values = fields[9 + sample_i].split(':')
            if len(values) < len(fmt):
                values.extend(['.'] * (len(fmt) - len(values)))
            gt = values[gt_idx]
            alleles = gt.replace('|', '/').split('/')
            gt_missing = any(a == '.' for a in alleles)
            if not gt_missing and len(alleles) != ploidy:
                raise SystemExit('GT ploidy mismatch at {}:{} sample {}: {}'.format(chrom, pos, sample, gt))
            if any(a not in ('.', '0', '1') for a in alleles):
                raise SystemExit('non-biallelic GT reached mask stage at {}:{} sample {}: {}'.format(chrom, pos, sample, gt))
            dp = parse_int(values[dp_idx], 'DP', chrom, pos, sample)
            gq = parse_int(values[gq_idx], 'GQ', chrom, pos, sample)
            fail_dp = dp is None or dp < 5
            fail_gq = gq is None or gq < 20
            if gt_missing:
                originally_missing += 1
                values[gt_idx] = canonical_missing(ploidy)
                miss_site += 1
            elif fail_dp or fail_gq:
                masked += 1
                low_dp += int(fail_dp)
                low_gq += int(fail_gq)
                missing_metric += int(dp is None or gq is None)
                values[gt_idx] = canonical_missing(ploidy)
                miss_site += 1
            else:
                # A valid male-X genotype must remain haploid; other valid calls are diploid.
                if len(alleles) != ploidy:
                    raise SystemExit('valid GT ploidy mismatch at {}:{} sample {}'.format(chrom, pos, sample))
            fields[9 + sample_i] = ':'.join(values)
        total += 1
        chrom_counts[chrom] += 1
        fmissing = miss_site / float(len(sexes))
        fmissing_sum += fmissing
        info_parts.append('F_MISSING={:.12g}'.format(fmissing))
        fields[7] = ';'.join(info_parts) if info_parts else '.'
        sys.stdout.write('\t'.join(fields) + '\n')

    if sample_names is None:
        raise SystemExit('VCF header was not found')
    sys.stderr.write('mask_summary records={} samples={} masked_low_or_missing_quality={} lowDP={} lowGQ={} missing_metric={} originally_missing_samples={} mean_F_MISSING={} chrom_counts={}\n'.format(
        total, len(sexes), masked, low_dp, low_gq, missing_metric, originally_missing,
        fmissing_sum / total if total else 0.0, dict(chrom_counts)))


if __name__ == '__main__':
    main()
