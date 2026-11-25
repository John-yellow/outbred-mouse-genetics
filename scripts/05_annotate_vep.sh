#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
source "$SCRIPT_DIR/lib/config.sh"
require_tool "$BCFTOOLS"
require_tool "$VEP_BIN"
ref=$(repo_path "$REFERENCE_FASTA")
cache=$(repo_path "$VEP_CACHE_DIR")
filter_dir=$(repo_path "$FILTER_DIR")
out=$(repo_path "$ANNOTATION_DIR")
mkdir -p "$out" "$out/logs"
candidate="$filter_dir/candidate_variants.bcf"
main="$filter_dir/main_variants.bcf"
for f in "$candidate" "$main"; do [[ -s "$f" ]] || { echo "Missing filtered callset $f" >&2; exit 2; }; done
candidate_vcf="$out/candidate_input.vcf.gz"
vep_site="$out/candidate_vep115_sites.vcf.gz"
candidate_out="$out/candidate_variants_vep115_GRCm39.bcf"
main_out="$out/main_variants_vep115_GRCm39.bcf"
for f in "$candidate_vcf" "$vep_site" "$candidate_out" "$main_out"; do
  [[ ! -e "$f" ]] || { echo "Refusing to overwrite $f" >&2; exit 2; }
done
"$BCFTOOLS" view -G -Oz -o "$candidate_vcf" "$candidate"
"$BCFTOOLS" index -c -f "$candidate_vcf"
"$VEP_BIN" --species mus_musculus --assembly GRCm39 --input_file "$candidate_vcf" \
  --format vcf --output_file "$vep_site" --compress_output bgzip --vcf --cache \
  --cache_version "$VEP_CACHE_VERSION" --dir_cache "$cache" --fasta "$ref" --offline \
  --everything --check_ref --check_existing --no_stats --fork "$VEP_FORKS" \
  --buffer_size "$VEP_BUFFER_SIZE" --warning_file "$out/logs/vep_warnings.log" \
  > "$out/logs/vep.stdout.log" 2> "$out/logs/vep.stderr.log"
"$BCFTOOLS" index -c -f "$vep_site"
for item in "candidate:$candidate:$candidate_out" "main:$main:$main_out"; do
  IFS=: read -r label src dest <<< "$item"
  "$BCFTOOLS" annotate --no-version --collapse none --threads 8 -a "$vep_site" \
    -c CHROM,POS,REF,ALT,INFO/CSQ -Ob -o "$dest" "$src"
  "$BCFTOOLS" view -h "$dest" >/dev/null
  "$BCFTOOLS" index -f "$dest"
  n_rec=$("$BCFTOOLS" index -n "$dest")
  [[ "$n_rec" =~ ^[0-9]+$ ]] || { echo "Invalid record count for $dest: $n_rec" >&2; exit 2; }
  "$BCFTOOLS" query -l "$src" > "$out/${label}.source.samples.txt"
  "$BCFTOOLS" query -l "$dest" > "$out/${label}.annotated.samples.txt"
  cmp "$out/${label}.source.samples.txt" "$out/${label}.annotated.samples.txt"
  header=$("$BCFTOOLS" view -h "$dest")
  [[ "$header" == *"##INFO=<ID=CSQ,"* ]] || { echo "Missing CSQ header in $dest" >&2; exit 2; }
done
