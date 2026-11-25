#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

require_cmd bcftools plink admixture awk sort

ANNOTATED_VCF="$(pipeline_path "${ANNOTATED_VCF:-data/vcf/obm2026_snps_annotated.vcf.gz}")"
RENAMED_VCF="$(pipeline_path "${RENAMED_VCF:-results/genetics/obm2026_snps_renamedID.vcf.gz}")"
GENETICS_DIR="$(pipeline_path "${GENETICS_DIR:-results/genetics}")"
GROUPS_FILE="$(pipeline_path "${SAMPLE_GROUPS:-config/sample_groups.tsv}")"
ADMIXTURE_K="${ADMIXTURE_K:-3}"
ADMIXTURE_CV="${ADMIXTURE_CV:-10}"
ADMIXTURE_THREADS="${ADMIXTURE_THREADS:-4}"

[[ -f "$ANNOTATED_VCF" ]] || { echo "ERROR: annotated VCF not found: $ANNOTATED_VCF" >&2; exit 1; }
[[ -f "$GROUPS_FILE" ]] || { echo "ERROR: sample-group file not found: $GROUPS_FILE" >&2; exit 1; }
mkdir -p "$GENETICS_DIR"

log "assigning variant IDs as CHROM:POS"
bcftools annotate --set-id '%CHROM:%POS' "$ANNOTATED_VCF" -Oz -o "$RENAMED_VCF"
bcftools index -t "$RENAMED_VCF"

BASE="$GENETICS_DIR/obm2026_snps_renamedID"
plink --vcf "$RENAMED_VCF" --make-bed --allow-extra-chr --out "$BASE"

# PLINK VCF imports do not reliably encode the biological population in FID.
# Update FID from the sample table so PCA/FST/ADMIXTURE downstream files carry
# explicit population labels. The first column in the table is the VCF/PLINK IID.
UPDATES="$GENETICS_DIR/sample_id_updates.txt"
MISSING="$GENETICS_DIR/sample_ids_missing_from_group_table.txt"
: > "$MISSING"
awk -v missing="$MISSING" '
BEGIN { FS = "[[:space:]]+"; OFS = "\t" }
FNR == NR {
    if ($0 !~ /^[[:space:]]*#/ && $1 != "sample_id" && NF >= 2) population[$1] = $2
    next
}
{
    if (!($2 in population)) print $2 > missing
    else print $1, $2, population[$2], $2
}
' "$GROUPS_FILE" "$BASE.fam" > "$UPDATES"

if [[ -s "$MISSING" ]]; then
    echo "ERROR: these PLINK sample IDs are absent from $GROUPS_FILE:" >&2
    sed -n '1,20p' "$MISSING" >&2
    exit 1
fi

GROUPED_BASE="$GENETICS_DIR/obm2026_snps_renamedID_grouped"
plink --bfile "$BASE" --update-ids "$UPDATES" --make-bed --allow-extra-chr --out "$GROUPED_BASE"
ANALYSIS_BASE="$GROUPED_BASE"

LD_BASE="$GENETICS_DIR/obm_ld_pruned"
PRUNED_BASE="$GENETICS_DIR/obm_pruned_dataset"
PCA_BASE="$GENETICS_DIR/obm_pca"

plink --bfile "$ANALYSIS_BASE" \
    --indep-pairwise 50 5 0.2 \
    --allow-extra-chr \
    --out "$LD_BASE"

plink --bfile "$ANALYSIS_BASE" \
    --extract "$LD_BASE.prune.in" \
    --make-bed \
    --allow-extra-chr \
    --out "$PRUNED_BASE"

plink --bfile "$PRUNED_BASE" \
    --pca \
    --allow-extra-chr \
    --out "$PCA_BASE"

CLUSTER_FILE="$GENETICS_DIR/population_clusters.txt"
awk '{print $1, $2, $1}' "$PRUNED_BASE.fam" > "$CLUSTER_FILE"
mapfile -t POPS < <(awk '{print $1}' "$PRUNED_BASE.fam" | sort -u)

FST_SUMMARY="$GENETICS_DIR/fst_summary.tsv"
printf 'Pop1\tPop2\tMean_Fst\tWeighted_Fst\n' > "$FST_SUMMARY"
for ((i = 0; i < ${#POPS[@]}; i++)); do
    for ((j = i + 1; j < ${#POPS[@]}; j++)); do
        P1="${POPS[i]}"
        P2="${POPS[j]}"
        FST_PREFIX="$GENETICS_DIR/fst_${P1}_vs_${P2}"
        log "FST: $P1 vs $P2"
        plink --bfile "$PRUNED_BASE" \
            --within "$CLUSTER_FILE" \
            --keep-cluster-names "$P1" "$P2" \
            --fst \
            --allow-extra-chr \
            --out "$FST_PREFIX"
        MEAN_FST="$(awk '/Mean Fst estimate:/ {print $4; exit}' "$FST_PREFIX.log")"
        WEIGHTED_FST="$(awk '/Weighted Fst estimate:/ {print $4; exit}' "$FST_PREFIX.log")"
        printf '%s\t%s\t%s\t%s\n' "$P1" "$P2" "$MEAN_FST" "$WEIGHTED_FST" >> "$FST_SUMMARY"
    done
done

(cd "$GENETICS_DIR" && admixture \
    --cv="$ADMIXTURE_CV" \
    --thread "$ADMIXTURE_THREADS" \
    "$(basename "$PRUNED_BASE").bed" "$ADMIXTURE_K")

log "stage 6 complete: $GENETICS_DIR"
