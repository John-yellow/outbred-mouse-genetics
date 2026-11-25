#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
source "$SCRIPT_DIR/lib/config.sh"
require_tool "$FAST_DUP"
manifest=$(repo_path "$SAMPLE_MANIFEST")
aligned=$(repo_path "$ALIGN_DIR")
marked=$(repo_path "$MARKED_DIR")
mkdir -p "$marked" "$marked/metrics"
header=1
while IFS=$'\t' read -r sample group sex r1 r2 extra; do
  if (( header )); then header=0; continue; fi
  [[ -z "${sample:-}" ]] && continue
  input="$aligned/${sample}.sorted.bam"
  output="$marked/${sample}.marked.bam"
  [[ -s "$input" ]] || { echo "Missing coordinate-sorted BAM: $input" >&2; exit 2; }
  [[ ! -e "$output" ]] || { echo "Refusing to overwrite $output" >&2; exit 2; }
  "$FAST_DUP" --input "$input" --metrics "$marked/metrics/${sample}.txt" \
    --output "$output" --num-threads "$DUP_THREADS" \
    --duplicate-scoring-strategy SUM_OF_BASE_QUALITIES \
    --read-name-regex null --create-index --index-format BAI
done < "$manifest"
