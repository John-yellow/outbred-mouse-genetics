#!/usr/bin/env bash

set -euo pipefail
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

require_cmd samtools bcftools parallel find awk sort

REF="$(pipeline_path "${REFERENCE_FASTA:-reference/Mus_musculus.GRCm39.dna.primary_assembly.fa}")"
BAM_DIR="$(pipeline_path "${BAM_DIR:-data/bam}")"
OUTPUT_DIR="$(pipeline_path "${VCF_DIR:-data/vcf}")"
WORK_DIR="$(pipeline_path "${WORK_DIR:-work}")"
LIST_FILE="$WORK_DIR/bams.list"
REGIONS_FILE="$WORK_DIR/regions.txt"
CONCAT_LIST="$WORK_DIR/vcf_concat_list.txt"
FINAL_OUT="$(pipeline_path "${RAW_VCF:-data/vcf/obm.raw.vcf.gz}")"

TARGET_CHROMS="${TARGET_CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y MT}"
BCFTOOLS_THREADS="${BCFTOOLS_THREADS:-2}"
CALL_JOBS="${CALL_JOBS:-4}"
REGION_SIZE="${REGION_SIZE:-10000000}"

[[ -f "$REF" ]] || { echo "ERROR: reference FASTA not found: $REF" >&2; exit 1; }
mkdir -p "$OUTPUT_DIR" "$WORK_DIR"
find "$BAM_DIR" -maxdepth 1 -type f -name '*.sorted.bam' | sort > "$LIST_FILE"

FIRST_BAM="$(head -n 1 "$LIST_FILE" || true)"
if [[ -z "$FIRST_BAM" || ! -f "$FIRST_BAM" ]]; then
    echo "ERROR: no sorted BAM files found under $BAM_DIR" >&2
    exit 1
fi

log "generating ${REGION_SIZE}-bp genomic windows"
awk -v size="$REGION_SIZE" -v target="$TARGET_CHROMS" '
BEGIN {
    n = split(target, wanted, /[[:space:]]+/)
    for (i = 1; i <= n; i++) keep[wanted[i]] = 1
}
$1 == "@SQ" {
    chrom = ""
    chrom_len = 0
    for (i = 2; i <= NF; i++) {
        if ($i ~ /^SN:/) chrom = substr($i, 4)
        if ($i ~ /^LN:/) chrom_len = substr($i, 4) + 0
    }
    if (keep[chrom] && chrom_len > 0) {
        for (start = 1; start <= chrom_len; start += size) {
            end = start + size - 1
            if (end > chrom_len) end = chrom_len
            print chrom ":" start "-" end
        }
    }
}
' < <(samtools view -H "$FIRST_BAM") > "$REGIONS_FILE"

[[ -s "$REGIONS_FILE" ]] || { echo "ERROR: no genomic regions were generated" >&2; exit 1; }

call_region() {
    set -o pipefail
    local region="$1"
    local safe_region="${region//:/_}"
    local output="$OUTPUT_DIR/obm.${safe_region}.vcf.gz"

    log "calling region: $region"
    bcftools mpileup \
        --threads "$BCFTOOLS_THREADS" \
        -g 10 \
        -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,INFO/AD \
        -E -Q 0 -pm 3 -F 0.25 -d 500 \
        -f "$REF" \
        -r "$region" \
        -b "$LIST_FILE" \
        -Ou | \
    bcftools call -mAv -f GQ,GP -p 0.99 -Oz -o "$output"
    bcftools index -t "$output"
}

export -f call_region log timestamp
export REF LIST_FILE OUTPUT_DIR BCFTOOLS_THREADS

log "running parallel SNP calling with $CALL_JOBS jobs"
parallel --will-cite --halt soon,fail=1 -j "$CALL_JOBS" call_region :::: "$REGIONS_FILE"

log "concatenating interval VCFs"
: > "$CONCAT_LIST"
while IFS= read -r region; do
    safe_region="${region//:/_}"
    printf '%s\n' "$OUTPUT_DIR/obm.${safe_region}.vcf.gz" >> "$CONCAT_LIST"
done < "$REGIONS_FILE"

bcftools concat -f "$CONCAT_LIST" -a -Oz -o "$FINAL_OUT"
bcftools index -t "$FINAL_OUT"

while IFS= read -r interval_vcf; do
    rm -f -- "$interval_vcf" "${interval_vcf}.tbi" "${interval_vcf}.csi"
done < "$CONCAT_LIST"
rm -f -- "$CONCAT_LIST" "$REGIONS_FILE"

log "stage 3 complete: $FINAL_OUT"
