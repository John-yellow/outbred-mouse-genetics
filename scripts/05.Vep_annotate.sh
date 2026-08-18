#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

require_cmd vep

INPUT_VCF="$(pipeline_path "${FILTERED_VCF:-data/vcf/obm_snps.vcf.gz}")"
OUTPUT_VCF="$(pipeline_path "${ANNOTATED_VCF:-data/vcf/obm2026_snps_annotated.vcf.gz}")"
VEP_CACHE_DIR="$(pipeline_path "${VEP_CACHE_DIR:-reference/vep_cache}")"
VEP_FASTA="$(pipeline_path "${VEP_FASTA:-reference/Mus_musculus.GRCm39.dna.primary_assembly.fa}")"
WORK_DIR="$(pipeline_path "${WORK_DIR:-work}")"
VEP_FORKS="${VEP_FORKS:-32}"

[[ -f "$INPUT_VCF" ]] || { echo "ERROR: input VCF not found: $INPUT_VCF" >&2; exit 1; }
[[ -f "$VEP_FASTA" ]] || { echo "ERROR: VEP FASTA not found: $VEP_FASTA" >&2; exit 1; }
[[ -d "$VEP_CACHE_DIR" ]] || { echo "ERROR: VEP cache directory not found: $VEP_CACHE_DIR" >&2; exit 1; }
mkdir -p "$(dirname "$OUTPUT_VCF")" "$WORK_DIR"

vep \
    --species mus_musculus \
    --assembly GRCm39 \
    --input_file "$INPUT_VCF" \
    --format vcf \
    --output_file "$OUTPUT_VCF" \
    --compress_output bgzip \
    --vcf \
    --cache \
    --cache_version 115 \
    --dir_cache "$VEP_CACHE_DIR" \
    --fasta "$VEP_FASTA" \
    --offline \
    --database 0 \
    --everything \
    --fork "$VEP_FORKS" \
    --warning_file "$WORK_DIR/vep.warnings.log" \
    --force_overwrite

if [[ -f "$OUTPUT_VCF" ]] && command -v bcftools >/dev/null 2>&1; then
    bcftools index -t "$OUTPUT_VCF" 2>/dev/null || true
fi
log "stage 5 complete: $OUTPUT_VCF"
