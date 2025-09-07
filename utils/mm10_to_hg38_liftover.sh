 #!/bin/bash

# LiftOver Script: Convert genomic coordinates from mm10 to hg38
# Uses UCSC liftOver tool to map mouse genome coordinates to human genome
# Version: 1.0 (Compatible with liftOver v377+)

set -euo pipefail

# Function to display usage
usage() {
    cat << EOF
Usage: $0 [input_bed] [chain_file] [output_success] [output_failed] [min_match]

Convert mm10 genomic coordinates to hg38 using UCSC liftOver.

Arguments (all optional, defaults shown):
  input_bed      Input BED file (results/mm10_merged_peaks_030925.bed)
  chain_file     Chain file (mm10ToHg38.over.chain.gz)
  output_success Output for successful conversions (results/liftover_mm10_to_hg38_peaks.bed)
  output_failed  Output for failed conversions (results/liftover_mm10_unsuccessful.bed)
  min_match      Minimum match ratio (0.1)

Options:
  -h, --help     Show this help message
  -v, --version  Show version information
EOF
}

# Handle help and version flags
case "${1:-}" in
    -h|--help) usage; exit 0 ;;
    -v|--version) echo "mm10_to_hg38_liftover v1.0"; exit 0 ;;
esac

# Default parameters
INPUT_BED="${1:-results/mm10_merged_peaks_030925.bed}"
CHAIN_FILE="${2:-mm10ToHg38.over.chain.gz}"
OUTPUT_SUCCESSFUL="${3:-results/liftover_mm10_to_hg38_peaks.bed}"
OUTPUT_FAILED="${4:-results/liftover_mm10_unsuccessful.bed}"
MIN_MATCH="${5:-0.1}"

# Validate required files
[[ -f "$INPUT_BED" ]] || { echo "Error: Input BED file not found: $INPUT_BED" >&2; exit 1; }
[[ -f "$CHAIN_FILE" ]] || { echo "Error: Chain file not found: $CHAIN_FILE" >&2; exit 1; }

# Create output directories if needed
mkdir -p "$(dirname "$OUTPUT_SUCCESSFUL")" "$(dirname "$OUTPUT_FAILED")"

# Execute liftOver
echo "Converting $INPUT_BED (mm10 → hg38)..."
liftOver -multiple -minMatch="$MIN_MATCH" "$INPUT_BED" "$CHAIN_FILE" "$OUTPUT_SUCCESSFUL" "$OUTPUT_FAILED"

# Report results
if [[ -f "$OUTPUT_SUCCESSFUL" ]] && [[ -f "$OUTPUT_FAILED" ]]; then
    successful=$(wc -l < "$OUTPUT_SUCCESSFUL")
    failed=$(wc -l < "$OUTPUT_FAILED")
    total=$((successful + failed))
    echo "Results: $successful/$total regions converted successfully ($(( successful * 100 / total ))%)"
    echo "Successful: $OUTPUT_SUCCESSFUL"
    echo "Failed: $OUTPUT_FAILED"
else
    echo "Error: liftOver execution failed" >&2
    exit 1
fi