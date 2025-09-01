#!/bin/bash

# Script: merge_bed_peaks.sh
# Description: Merges BED files, generating a consolidated BED file and a statistics log.
# Usage: bash merge_bed_peaks.sh <input_bed_file1> <input_bed_file2> ... <input_bed_filen-1> <output_name>

# $1-$(n-1) list of bed files to merge
# $n - output file name (bed extension will be added automatically)

# Check if at least one file is provided as an argument
if [ "$#" -lt 1 ]; then
  echo "No files provided."
  exit 1
fi

output="${!#}"
temp_dir='tmp'

mkdir -p "$temp_dir"

log_file="$output.log"
total_lines=0
total_nucleotides=0

# Clear the log file if it already exists
> "$log_file"

# Output file
output_file="tmp/combined.bed"


# Concatenate BED files
for file in "$@"; do
  if [ -f "$file" ]; then
    # Count the number of lines
    cat "$file" >> "$output_file"
    line_count=$(wc -l < "$file")
    total_lines=$((total_lines + line_count))

    # Calculate the sum of nucleotides
    nucleotides_sum=$(awk -F'\t' 'BEGIN{diff_sum = 0} {diff = $3 - $2; diff_sum += diff} END{print diff_sum}' "$file")

    # Write the file name, line count, and nucleotide sum to the log file
    echo -e "$file\t$line_count\t$nucleotides_sum" >> "$log_file"
  fi
done

# sort all the peaks
sorted_combined="tmp/sort_combined.bed"
sort -k1,1 -k2,2n "$output_file" > "$sorted_combined"
rm -f "$output_file"

# merge 
bedtools merge -i "$sorted_combined" -c 4 -o distinct > "$output.bed"

after_merge_line_count=$(wc -l < "$output.bed")
gap_peaks=$((total_lines - after_merge_line_count))
sum_before=$(awk -F'\t' 'BEGIN{diff_sum = 0} {diff = $3 - $2; diff_sum += diff} END{print diff_sum}' "$output.bed")

# sum the number of bp in the merged file
sum_after=$(awk -F'\t' 'BEGIN{diff_sum = 0} {diff = $3 - $2; diff_sum += diff} END{print diff_sum}' "$output.bed")

gap_bp=$((sum_before - sum_after))

echo -e "# peaks before merging - $total_lines\n# bp before - $sum_before \n# peaks after merging - $after_merge_line_count\n# bp after - $sum_after\nGap of peaks - $gap_peaks\nGap of bp -  $gap_bp" >> "$log_file"

rm -f "$sorted_combined"
rmdir "$temp_dir"

cat "$log_file"
echo Done!
