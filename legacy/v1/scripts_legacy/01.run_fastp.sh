#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

require_cmd fastp xargs find awk

RAW_DIR="$(pipeline_path "${RAW_DIR:-data/raw}")"
CLEAN_DIR="$(pipeline_path "${CLEAN_DIR:-data/clean}")"
MAP_FILE="$(pipeline_path "${SAMPLE_MAP:-config/sample_map.tsv}")"
FASTP_THREADS="${FASTP_THREADS:-3}"
FASTP_JOBS="${FASTP_JOBS:-4}"

mkdir -p "$CLEAN_DIR"

if [[ ! -f "$MAP_FILE" ]]; then
    log "sample map not found; generating a two-column map from raw directories"
    mkdir -p "$(dirname "$MAP_FILE")"
    find "$RAW_DIR" -mindepth 1 -maxdepth 1 -type d -printf '%f\t%f\n' | sort > "$MAP_FILE"
fi

if [[ ! -d "$RAW_DIR" ]]; then
    echo "ERROR: raw-read directory does not exist: $RAW_DIR" >&2
    exit 1
fi

process_sample() {
    local raw_name="$1"
    local clean_name="$2"
    local sample_dir="$RAW_DIR/$raw_name"
    local target_dir="$CLEAN_DIR/$clean_name"
    local r1_file
    local r2_file

    if [[ ! -d "$sample_dir" ]]; then
        echo "ERROR: sample directory does not exist: $sample_dir" >&2
        return 1
    fi

    if [[ -s "$target_dir/${clean_name}_R1_clean.fq.gz" && \
          -s "$target_dir/${clean_name}_R2_clean.fq.gz" ]]; then
        log "skip: $clean_name (clean reads already exist)"
        return 0
    fi

    r1_file="$(find "$sample_dir" -maxdepth 1 -type f \( \
        -name '*R1*.fq.gz' -o -name '*R1*.fastq.gz' -o \
        -name '*_1.fq.gz' -o -name '*_1.fastq.gz' -o \
        -name '*R1*.fq' -o -name '*R1*.fastq' -o \
        -name '*_1.fq' -o -name '*_1.fastq' \) | sort | head -n 1)"
    r2_file="$(find "$sample_dir" -maxdepth 1 -type f \( \
        -name '*R2*.fq.gz' -o -name '*R2*.fastq.gz' -o \
        -name '*_2.fq.gz' -o -name '*_2.fastq.gz' -o \
        -name '*R2*.fq' -o -name '*R2*.fastq' -o \
        -name '*_2.fq' -o -name '*_2.fastq' \) | sort | head -n 1)"

    if [[ -z "$r1_file" || -z "$r2_file" ]]; then
        echo "ERROR: paired R1/R2 reads not found for $raw_name in $sample_dir" >&2
        return 1
    fi

    mkdir -p "$target_dir"
    log "fastp: $raw_name -> $clean_name"
    fastp \
        --in1 "$r1_file" \
        --in2 "$r2_file" \
        --out1 "$target_dir/${clean_name}_R1_clean.fq.gz" \
        --out2 "$target_dir/${clean_name}_R2_clean.fq.gz" \
        --html "$target_dir/${clean_name}_report.html" \
        --json "$target_dir/${clean_name}_report.json" \
        --thread "$FASTP_THREADS" \
        --compression 1
    log "finished: $clean_name"
}

export -f process_sample
export RAW_DIR CLEAN_DIR FASTP_THREADS

if ! awk 'BEGIN {FS="[[:space:]]+"} $0 !~ /^[[:space:]]*#/ && NF >= 2 {print $1, $2}' "$MAP_FILE" | \
    xargs -r -n 2 -P "$FASTP_JOBS" bash -c 'process_sample "$1" "$2"' _; then
    echo "ERROR: at least one fastp job failed" >&2
    exit 1
fi

log "stage 1 complete"
