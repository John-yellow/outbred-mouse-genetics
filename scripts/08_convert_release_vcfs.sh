#!/usr/bin/env bash
# Convert annotated BCF callsets to BGZF-compressed VCF.gz with CSI indexes and validate output integrity.
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
source "$SCRIPT_DIR/lib/config.sh"
require_tool "$BCFTOOLS"
annotation=$(repo_path "$ANNOTATION_DIR")
staging=$(repo_path "${RELEASE_DIR:-release_staging}")
threads="${RELEASE_THREADS:-8}"
sex_map="${SEX_MAP:-}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads) threads="$2"; shift 2 ;;
    --sex-map) sex_map="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done
mkdir -p "$staging"
for set in candidate main; do
  if [[ "$set" == candidate ]]; then
    source="$annotation/candidate_variants_vep115_GRCm39.bcf"
    target="$staging/candidate_variants_vep115_GRCm39.vcf.gz"
  else
    source="$annotation/main_variants_vep115_GRCm39.bcf"
    target="$staging/main_variants_vep115_GRCm39.vcf.gz"
  fi
  [[ -s "$source" ]] || { echo "Missing annotated source BCF: $source" >&2; exit 2; }
  [[ ! -e "$target" && ! -e "$target.csi" ]] || { echo "Refusing to overwrite $target" >&2; exit 2; }
  "$BCFTOOLS" view --threads "$threads" --no-version -Oz -o "$target" "$source"
  "$BCFTOOLS" index --csi -f "$target"
  validate_args=(
    --source-bcf "$source"
    --vcf-gz "$target"
    --bcftools "$BCFTOOLS"
    --report "$staging/${set}_vcf_validation.json"
  )
  if [[ -n "$sex_map" ]]; then
    validate_args+=(--sex-map "$(repo_path "$sex_map")")
  fi
  python3 "$SCRIPT_DIR/validate_release_vcf.py" "${validate_args[@]}"
done
