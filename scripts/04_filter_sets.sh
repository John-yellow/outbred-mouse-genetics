#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
source "$SCRIPT_DIR/lib/config.sh"
require_tool "$BCFTOOLS"
ref=$(repo_path "$REFERENCE_FASTA")
manifest=$(repo_path "$SAMPLE_MANIFEST")
call=$(repo_path "$CALL_DIR")
out=$(repo_path "$FILTER_DIR")
mkdir -p "$out" "$call"
raw="$call/joint_call.raw.bcf"
sex_map="$call/sex.tsv"
[[ -s "$raw" ]] || { echo "Missing $raw; run scripts/03_call_windows.py first" >&2; exit 2; }
awk -F '\t' 'NR>1 {print $1 "\t" $3}' "$manifest" > "$sex_map"
normalized="$out/normalized.bcf"
masked="$out/masked.bcf"
trimmed="$out/trimmed.bcf"
with_tags="$out/genotype_qc.bcf"
candidate="$out/candidate_variants.bcf"
main="$out/main_variants.bcf"
for f in "$normalized" "$masked" "$trimmed" "$with_tags" "$candidate" "$main"; do
  [[ ! -e "$f" ]] || { echo "Refusing to overwrite $f" >&2; exit 2; }
done
"$BCFTOOLS" norm -f "$ref" -c e -m -any -Ob -o "$normalized" "$raw"
"$BCFTOOLS" view -Ov "$normalized" \
  | python3 "$SCRIPT_DIR/mask_genotypes.py" --sex-map "$sex_map" \
  | "$BCFTOOLS" view -Ob -o "$masked"
"$BCFTOOLS" view -a -Ob -o "$trimmed" "$masked"
"$BCFTOOLS" +fill-tags "$trimmed" -Ob -o "$with_tags" -- -t AC,AN,AF,MAF,NS
candidate_expr='QUAL>=30 && INFO/AC>0 && INFO/AC<INFO/AN'
"$BCFTOOLS" view -m 2 -M 2 -v snps -i "$candidate_expr" -Ob -o "$candidate" "$with_tags"
main_expr="$candidate_expr && INFO/F_MISSING<=0.05"
"$BCFTOOLS" view -m 2 -M 2 -v snps -i "$main_expr" -Ob -o "$main" "$with_tags"
expected_samples=$(wc -l < "$sex_map")
for f in "$candidate" "$main"; do
  "$BCFTOOLS" view -h "$f" >/dev/null
  "$BCFTOOLS" index -f "$f"
  n_rec=$("$BCFTOOLS" index -n "$f")
  [[ "$n_rec" =~ ^[0-9]+$ ]] || { echo "Invalid record count for $f: $n_rec" >&2; exit 2; }
  n_samp=$("$BCFTOOLS" query -l "$f" | wc -l)
  [[ "$n_samp" -eq "$expected_samples" ]] || { echo "Sample count mismatch in $f: $n_samp != $expected_samples" >&2; exit 2; }
  echo "$f: records=$n_rec samples=$n_samp"
done
