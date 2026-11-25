#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

require_cmd plink awk

RENAMED_VCF="$(pipeline_path "${RENAMED_VCF:-results/genetics/obm2026_snps_renamedID.vcf.gz}")"
GENETICS_DIR="$(pipeline_path "${GENETICS_DIR:-results/genetics}")"
GROUPS_FILE="$(pipeline_path "${SAMPLE_GROUPS:-config/sample_groups.tsv}")"
POPULATION_DIR="$(pipeline_path "${POPULATION_DIR:-results/genetics/populations}")"

[[ -f "$RENAMED_VCF" ]] || { echo "ERROR: renamed VCF not found: $RENAMED_VCF" >&2; exit 1; }
[[ -f "$GROUPS_FILE" ]] || { echo "ERROR: sample-group file not found: $GROUPS_FILE" >&2; exit 1; }
mkdir -p "$GENETICS_DIR" "$POPULATION_DIR"

BASE="$GENETICS_DIR/obm2026_snps_renamedID"
if [[ ! -s "$BASE.bed" || ! -s "$BASE.bim" || ! -s "$BASE.fam" ]]; then
    plink --vcf "$RENAMED_VCF" --make-bed --allow-extra-chr --out "$BASE"
fi

for population in ICR1 ICR2 KM; do
    KEEP_FILE="$POPULATION_DIR/${population}.keep"
    OUTPUT_BASE="$POPULATION_DIR/$population"

    awk -v population="$population" '
    BEGIN { FS = "[[:space:]]+"; OFS = "\t" }
    FNR == NR {
        if ($0 !~ /^[[:space:]]*#/ && $1 != "sample_id" && NF >= 2) group[$1] = $2
        next
    }
    group[$2] == population { print $1, $2 }
    ' "$GROUPS_FILE" "$BASE.fam" > "$KEEP_FILE"

    if [[ ! -s "$KEEP_FILE" ]]; then
        echo "ERROR: no PLINK samples assigned to population $population" >&2
        exit 1
    fi

    plink --bfile "$BASE" \
        --keep "$KEEP_FILE" \
        --maf 0.01 \
        --recode \
        --allow-extra-chr \
        --out "$OUTPUT_BASE"
done

log "stage 7 complete: $POPULATION_DIR"
