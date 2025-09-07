#!/bin/bash

# Script: compare_bed_files.sh
# Description: Compares two BED files and generates various intersection outputs
# Usage: bash compare_bed_files.sh <merged_bed> <reference_bed_with_ids> <output_dir>

# Input validation
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <merged_bed> <reference_bed_with_ids> <output_dir>"
    echo "Example: $0 merged_mice_mm10_atac.bed mouse_mm10_peaks_withID_1406.bed comparison_results"
    exit 1
fi

MERGED_BED="$1"
REFERENCE_BED="$2"
OUTPUT_DIR="$3"

# Check if input files exist
if [ ! -f "$MERGED_BED" ]; then
    echo "Error: Merged BED file '$MERGED_BED' not found."
    exit 1
fi

if [ ! -f "$REFERENCE_BED" ]; then
    echo "Error: Reference BED file '$REFERENCE_BED' not found."
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Output files
COMPLETE_OVERLAPS="$OUTPUT_DIR/complete_overlaps_with_ids.bed"
INTERSECTIONS="$OUTPUT_DIR/intersections_with_ids.bed"
MERGED_ONLY="$OUTPUT_DIR/merged_only_segments.bed"
ALL_COMBINED="$OUTPUT_DIR/all_segments_combined.bed"
STATS_FILE="$OUTPUT_DIR/comparison_statistics.txt"

echo "Starting BED file comparison..."
echo "Merged BED: $MERGED_BED"
echo "Reference BED: $REFERENCE_BED"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Clear output files
> "$COMPLETE_OVERLAPS"
> "$INTERSECTIONS"
> "$MERGED_ONLY"
> "$ALL_COMBINED"
> "$STATS_FILE"

echo "=== BED File Comparison Statistics ===" > "$STATS_FILE"
echo "Generated on: $(date)" >> "$STATS_FILE"
echo "Merged BED file: $MERGED_BED" >> "$STATS_FILE"
echo "Reference BED file: $REFERENCE_BED" >> "$STATS_FILE"
echo "" >> "$STATS_FILE"

# Count input lines
merged_count=$(wc -l < "$MERGED_BED")
reference_count=$(wc -l < "$REFERENCE_BED")

echo "Input file statistics:" >> "$STATS_FILE"
echo "  Merged BED peaks: $merged_count" >> "$STATS_FILE"
echo "  Reference BED peaks: $reference_count" >> "$STATS_FILE"
echo "" >> "$STATS_FILE"

# 1. Find complete overlaps (100% overlap in both directions)
echo "Finding complete overlaps..."
bedtools intersect -a "$REFERENCE_BED" -b "$MERGED_BED" -f 1.0 -r > "$COMPLETE_OVERLAPS"
complete_overlaps_count=$(wc -l < "$COMPLETE_OVERLAPS")

# 2. Find all intersections (any overlap) with IDs from reference
echo "Finding all intersections..."
bedtools intersect -a "$REFERENCE_BED" -b "$MERGED_BED" -wa > "$INTERSECTIONS"
intersections_count=$(wc -l < "$INTERSECTIONS")

# 3. Find segments only in merged file (not in reference)
echo "Finding merged-only segments..."
bedtools intersect -a "$MERGED_BED" -b "$REFERENCE_BED" -v > "$MERGED_ONLY"
merged_only_count=$(wc -l < "$MERGED_ONLY")

# 4. Create combined file with all segments
echo "Creating combined file..."
# First add all reference segments with their IDs
cat "$REFERENCE_BED" > "$ALL_COMBINED"
# Then add merged-only segments (those not overlapping with reference)
cat "$MERGED_ONLY" >> "$ALL_COMBINED"

# Sort the combined file
sort -k1,1 -k2,2n "$ALL_COMBINED" > "${ALL_COMBINED}.tmp"
mv "${ALL_COMBINED}.tmp" "$ALL_COMBINED"

all_combined_count=$(wc -l < "$ALL_COMBINED")

# Calculate statistics
echo "Output file statistics:" >> "$STATS_FILE"
echo "  Complete overlaps: $complete_overlaps_count" >> "$STATS_FILE"
echo "  All intersections: $intersections_count" >> "$STATS_FILE"
echo "  Merged-only segments: $merged_only_count" >> "$STATS_FILE"
echo "  All segments combined: $all_combined_count" >> "$STATS_FILE"
echo "" >> "$STATS_FILE"

# Calculate overlap percentages
if [ "$reference_count" -gt 0 ]; then
    complete_overlap_pct=$(awk "BEGIN {printf \"%.2f\", ($complete_overlaps_count/$reference_count)*100}")
    intersection_pct=$(awk "BEGIN {printf \"%.2f\", ($intersections_count/$reference_count)*100}")
else
    complete_overlap_pct="N/A"
    intersection_pct="N/A"
fi

if [ "$merged_count" -gt 0 ]; then
    merged_only_pct=$(awk "BEGIN {printf \"%.2f\", ($merged_only_count/$merged_count)*100}")
else
    merged_only_pct="N/A"
fi

echo "Overlap analysis:" >> "$STATS_FILE"
echo "  Complete overlap rate: $complete_overlap_pct% of reference peaks" >> "$STATS_FILE"
echo "  Intersection rate: $intersection_pct% of reference peaks" >> "$STATS_FILE"
echo "  Novel segments rate: $merged_only_pct% of merged peaks" >> "$STATS_FILE"
echo "" >> "$STATS_FILE"

# Calculate base pair statistics
echo "Calculating base pair coverage..."
reference_bp=$(awk -F'\t' '{sum += $3 - $2} END {print sum}' "$REFERENCE_BED")
merged_bp=$(awk -F'\t' '{sum += $3 - $2} END {print sum}' "$MERGED_BED")
complete_overlaps_bp=$(awk -F'\t' '{sum += $3 - $2} END {print sum}' "$COMPLETE_OVERLAPS")
intersections_bp=$(awk -F'\t' '{sum += $3 - $2} END {print sum}' "$INTERSECTIONS")
merged_only_bp=$(awk -F'\t' '{sum += $3 - $2} END {print sum}' "$MERGED_ONLY")

echo "Base pair coverage:" >> "$STATS_FILE"
echo "  Reference BED: $reference_bp bp" >> "$STATS_FILE"
echo "  Merged BED: $merged_bp bp" >> "$STATS_FILE"
echo "  Complete overlaps: $complete_overlaps_bp bp" >> "$STATS_FILE"
echo "  All intersections: $intersections_bp bp" >> "$STATS_FILE"
echo "  Merged-only segments: $merged_only_bp bp" >> "$STATS_FILE"
echo "" >> "$STATS_FILE"

# Generate file descriptions
echo "Output file descriptions:" >> "$STATS_FILE"
echo "  $COMPLETE_OVERLAPS - Peaks that completely overlap (100%) between both files" >> "$STATS_FILE"
echo "  $INTERSECTIONS - All peaks from reference that have any overlap with merged file" >> "$STATS_FILE"
echo "  $MERGED_ONLY - Peaks that exist only in merged file (novel segments)" >> "$STATS_FILE"
echo "  $ALL_COMBINED - Combined file with all segments from both sources" >> "$STATS_FILE"
echo "" >> "$STATS_FILE"

echo "Comparison completed successfully!"
echo ""
echo "=== Summary ==="
echo "Complete overlaps: $complete_overlaps_count peaks ($complete_overlap_pct% of reference)"
echo "All intersections: $intersections_count peaks ($intersection_pct% of reference)"
echo "Merged-only segments: $merged_only_count peaks ($merged_only_pct% of merged)"
echo "Total combined segments: $all_combined_count peaks"
echo ""
echo "Output files created in: $OUTPUT_DIR"
echo "See $STATS_FILE for detailed statistics"

