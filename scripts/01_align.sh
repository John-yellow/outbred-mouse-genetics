#!/usr/bin/env bash
# Default workflow input: provider-supplied quality-controlled paired-end reads.
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
source "$SCRIPT_DIR/lib/config.sh"
require_tool "$BWA_MEM2"
require_tool "$SAMTOOLS"
manifest=$(repo_path "$SAMPLE_MANIFEST")
ref=$(repo_path "$REFERENCE_FASTA")
out=$(repo_path "$ALIGN_DIR")
mkdir -p "$out" "$(dirname -- "$ref")"
[[ -s "$manifest" && -s "$ref" ]] || { echo "Manifest or reference is missing" >&2; exit 2; }
if [[ -n "${REFERENCE_SHA256:-}" ]]; then
  actual=$(sha256sum "$ref" | cut -d' ' -f1)
  [[ "$actual" == "$REFERENCE_SHA256" ]] || { echo "Reference SHA-256 mismatch" >&2; exit 2; }
fi
"$SAMTOOLS" faidx "$ref"
if [[ ! -s "${ref}.bwt.2bit.64" ]]; then
  "$BWA_MEM2" index "$ref"
fi
header=1
while IFS=$'\t' read -r sample group sex r1 r2 extra; do
  if (( header )); then header=0; continue; fi
  [[ -z "${sample:-}" ]] && continue
  [[ -z "${extra:-}" && "$sample" =~ ^[A-Za-z0-9_.-]+$ && "$sex" =~ ^(M|F)$ ]] || { echo "Invalid manifest row for $sample" >&2; exit 2; }
  r1=$(repo_path "$r1"); r2=$(repo_path "$r2")
  [[ -s "$r1" && -s "$r2" ]] || { echo "Missing paired reads for $sample" >&2; exit 2; }
  target="$out/${sample}.sorted.bam"
  partial="$target.partial.bam"
  [[ ! -e "$target" && ! -e "$partial" ]] || { echo "Refusing to overwrite $target" >&2; exit 2; }
  rg=$(printf '@RG\tID:%s\tSM:%s\tPL:ILLUMINA' "$sample" "$sample")
  "$BWA_MEM2" mem -t "$ALIGN_THREADS" -K 50000000 -R "$rg" "$ref" "$r1" "$r2" \
    | "$SAMTOOLS" sort -@ "$SORT_THREADS" -m "$SORT_MEMORY" -l "$SORT_LEVEL" -O BAM -o "$partial" -
  mv "$partial" "$target"
  "$SAMTOOLS" index -@ 4 "$target"
done < "$manifest"
