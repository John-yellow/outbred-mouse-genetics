#!/usr/bin/env bash

set -euo pipefail
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

require_cmd bwa samtools xargs find nproc

CLEAN_DIR="$(pipeline_path "${CLEAN_DIR:-data/clean}")"
BAM_DIR="$(pipeline_path "${BAM_DIR:-data/bam}")"
REF_GENOME="$(pipeline_path "${REFERENCE_FASTA:-reference/Mus_musculus.GRCm39.dna.primary_assembly.fa}")"
TARGET_CHROMS="${TARGET_CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y MT}"

BWA_THREADS="${BWA_THREADS:-8}"
VIEW_THREADS="${VIEW_THREADS:-1}"
SORT_THREADS="${SORT_THREADS:-2}"
MAPPING_JOBS="${MAPPING_JOBS:-1}"

[[ -f "$REF_GENOME" ]] || { echo "ERROR: reference FASTA not found: $REF_GENOME" >&2; exit 1; }
[[ -d "$CLEAN_DIR" ]] || { echo "ERROR: clean-read directory not found: $CLEAN_DIR" >&2; exit 1; }
mkdir -p "$BAM_DIR"

process_sample() {
    local sample_dir="$1"
    local r1_file
    local r2_file
    local sample_name
    local final_bam
    local bwa_log
    local chroms=()

    r1_file="$(find "$sample_dir" -maxdepth 1 -type f -name '*_R1_clean.fq.gz' | sort | head -n 1)"
    r2_file="$(find "$sample_dir" -maxdepth 1 -type f -name '*_R2_clean.fq.gz' | sort | head -n 1)"

    if [[ -z "$r1_file" || -z "$r2_file" ]]; then
        echo "ERROR: paired clean reads not found in $sample_dir" >&2
        return 1
    fi

    sample_name="$(basename "$r1_file" '_R1_clean.fq.gz')"
    final_bam="$BAM_DIR/${sample_name}.sorted.bam"
    bwa_log="$BAM_DIR/${sample_name}.bwa.log"
    read -r -a chroms <<< "$TARGET_CHROMS"

    if [[ -s "$final_bam" && -s "${final_bam}.bai" ]]; then
        if samtools quickcheck "$final_bam" >/dev/null 2>&1; then
            log "skip: $sample_name (valid BAM already exists)"
            return 0
        fi
        log "invalid BAM found; rebuilding: $sample_name"
        rm -f -- "$final_bam" "${final_bam}.bai"
    elif [[ -e "$final_bam" || -e "${final_bam}.bai" ]]; then
        log "incomplete BAM pair found; rebuilding: $sample_name"
        rm -f -- "$final_bam" "${final_bam}.bai"
    fi

    log "mapping: $sample_name"
    if ! bwa mem -t "$BWA_THREADS" \
        -R "@RG\tID:${sample_name}\tPL:ILLUMINA\tSM:${sample_name}" \
        "$REF_GENOME" "$r1_file" "$r2_file" 2>"$bwa_log" | \
        samtools view -b -@ "$VIEW_THREADS" - "${chroms[@]}" | \
        samtools sort -@ "$SORT_THREADS" -m 2G -o "$final_bam" -; then
        echo "ERROR: alignment failed for $sample_name; see $bwa_log" >&2
        rm -f -- "$final_bam" "${final_bam}.bai"
        return 1
    fi

    if [[ ! -s "$final_bam" ]] || ! samtools quickcheck "$final_bam" >/dev/null 2>&1; then
        echo "ERROR: invalid BAM produced for $sample_name; see $bwa_log" >&2
        rm -f -- "$final_bam" "${final_bam}.bai"
        return 1
    fi

    if ! samtools index -@ 4 "$final_bam"; then
        echo "ERROR: BAM indexing failed for $sample_name" >&2
        rm -f -- "$final_bam" "${final_bam}.bai"
        return 1
    fi
    log "finished: ${sample_name}.sorted.bam"
}

export -f process_sample
export REF_GENOME BAM_DIR TARGET_CHROMS BWA_THREADS VIEW_THREADS SORT_THREADS

sample_count="$(find "$CLEAN_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)"
if [[ "$sample_count" -eq 0 ]]; then
    echo "ERROR: no sample directories found under $CLEAN_DIR" >&2
    exit 1
fi

if ! find "$CLEAN_DIR" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z | \
    xargs -0 -r -n 1 -P "$MAPPING_JOBS" bash -c 'process_sample "$1"' _; then
    echo "ERROR: at least one mapping job failed" >&2
    exit 1
fi

log "stage 2 complete"
