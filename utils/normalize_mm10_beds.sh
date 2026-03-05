#!/usr/bin/env bash

set -euo pipefail

RAW_DIR="data/mm10_atac_seq/raw_files"
OUT_DIR="data/mm10_atac_seq"

# Ensure output directory exists
mkdir -p "$OUT_DIR"

# If no files match, the loop should not iterate
shopt -s nullglob

for bed in "$RAW_DIR"/*.bed; do
  base_name="$(basename "$bed" .bed)"
  out_path="$OUT_DIR/${base_name}_normalized.bed"
  # Keep only chr rows and the first 3 columns
  cut -f1-3 "$bed" | grep -E '^chr' > "$out_path"
  echo "Wrote: $out_path"
done

shopt -u nullglob


