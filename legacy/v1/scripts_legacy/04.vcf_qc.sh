#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

require_cmd bcftools

INPUT_VCF="$(pipeline_path "${RAW_VCF:-data/vcf/obm.raw.vcf.gz}")"
REF_FASTA="$(pipeline_path "${REFERENCE_FASTA:-reference/Mus_musculus.GRCm39.dna.primary_assembly.fa}")"
OUTPUT_VCF="$(pipeline_path "${FILTERED_VCF:-data/vcf/obm_snps.vcf.gz}")"
WORK_DIR="$(pipeline_path "${WORK_DIR:-work}/vcf_qc")"
TARGET_CHRS="${QC_CHROMS:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X}"

[[ -f "$INPUT_VCF" ]] || { echo "ERROR: input VCF not found: $INPUT_VCF" >&2; exit 1; }
[[ -f "$REF_FASTA" ]] || { echo "ERROR: reference FASTA not found: $REF_FASTA" >&2; exit 1; }
mkdir -p "$(dirname "$OUTPUT_VCF")" "$WORK_DIR"

STEP1="$WORK_DIR/step1.chr_selected.bcf"
STEP2="$WORK_DIR/step2.normed.bcf"
STEP3="$WORK_DIR/step3.recalculate.bcf"
STEP4="$WORK_DIR/step4.trim_alt.bcf"
STEP5="$WORK_DIR/step5.filtered.bcf"

log "selecting chromosomes: $TARGET_CHRS"
bcftools view -t "$TARGET_CHRS" -Ob -o "$STEP1" "$INPUT_VCF"
bcftools norm -m +any -f "$REF_FASTA" -Ob -o "$STEP2" "$STEP1"
bcftools plugin fill-tags "$STEP2" -Ob -o "$STEP3" -- -t AN,AC,MAF
bcftools view -a "$STEP3" -Ob -o "$STEP4"
bcftools view \
    -m2 -M2 -v snps \
    -i 'QUAL>=30 && MAF>0.01 && F_MISSING==0 && DP>=5 && GQ>=20' \
    -Ob -o "$STEP5" "$STEP4"
bcftools view -Oz -o "$OUTPUT_VCF" "$STEP5"
bcftools index -t "$OUTPUT_VCF"

record_count="$(bcftools index -n "$OUTPUT_VCF")"
log "stage 4 complete: $OUTPUT_VCF ($record_count records)"
