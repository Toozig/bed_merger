#!/bin/bash

# Script: bed_analysis_pipeline.sh
# Description: Comprehensive BED file analysis pipeline with GTF coding region extraction and Excel reporting
# Usage: bash bed_analysis_pipeline.sh <bed_files_directory> <gtf_file> <output_directory>

set -e  # Exit on any error

# Input validation
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <bed_files_directory> <gtf_file> <output_directory>"
    echo "Example: $0 mice_mm10_atac_seq_bed_Files/ mm10_genes.gtf analysis_output"
    exit 1
fi

BED_DIR="$1"
GTF_FILE="$2"
OUTPUT_DIR="$3"

# Get absolute paths
BED_DIR=$(readlink -f "$BED_DIR")
GTF_FILE=$(readlink -f "$GTF_FILE")
OUTPUT_DIR=$(readlink -f "$OUTPUT_DIR")

# Check if input files/directories exist
if [ ! -d "$BED_DIR" ]; then
    echo "Error: BED files directory '$BED_DIR' not found."
    exit 1
fi

if [ ! -f "$GTF_FILE" ]; then
    echo "Error: GTF file '$GTF_FILE' not found."
    exit 1
fi

# Create organized output directory structure
mkdir -p "$OUTPUT_DIR"

# Set up organized directory structure
RESULTS_DIR="$OUTPUT_DIR/results"
INPUT_COPY_DIR="$OUTPUT_DIR/input_files"
SCRIPTS_COPY_DIR="$OUTPUT_DIR/scripts"
TEMP_DIR="$OUTPUT_DIR/temp"

mkdir -p "$RESULTS_DIR" "$INPUT_COPY_DIR" "$SCRIPTS_COPY_DIR" "$TEMP_DIR"

# Set up file paths
CODING_BED="$RESULTS_DIR/coding_regions.bed"
MERGED_BED="$RESULTS_DIR/merged_bed_files.bed"
MERGED_LOG="$RESULTS_DIR/merged_bed_files.log"
EXCEL_OUTPUT="$RESULTS_DIR/bed_analysis_report.xlsx"
CONFIG_LOG="$OUTPUT_DIR/analysis_config.log"

# Record start time for timing analysis
start_time=$(date +%s)

echo "=== BED Analysis Pipeline ==="
echo "BED files directory: $BED_DIR"
echo "GTF file: $GTF_FILE"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Generate configuration log
echo "=== BED Analysis Pipeline Configuration Log ===" > "$CONFIG_LOG"
echo "Generated on: $(date)" >> "$CONFIG_LOG"
echo "Pipeline version: 1.0" >> "$CONFIG_LOG"
echo "User: $(whoami)" >> "$CONFIG_LOG"
echo "Host: $(hostname)" >> "$CONFIG_LOG"
echo "Working directory: $(pwd)" >> "$CONFIG_LOG"
echo "" >> "$CONFIG_LOG"

echo "=== Input Configuration ===" >> "$CONFIG_LOG"
echo "BED files directory: $BED_DIR" >> "$CONFIG_LOG"
echo "GTF file: $GTF_FILE" >> "$CONFIG_LOG"
echo "Output directory: $OUTPUT_DIR" >> "$CONFIG_LOG"
echo "" >> "$CONFIG_LOG"

echo "=== Software Versions ===" >> "$CONFIG_LOG"
echo "Bash: $BASH_VERSION" >> "$CONFIG_LOG"
python3 --version >> "$CONFIG_LOG" 2>&1
bedtools --version >> "$CONFIG_LOG" 2>&1
awk --version | head -1 >> "$CONFIG_LOG" 2>&1
sort --version | head -1 >> "$CONFIG_LOG" 2>&1

# Check Python packages
echo "" >> "$CONFIG_LOG"
echo "=== Python Package Versions ===" >> "$CONFIG_LOG"
python3 -c "import pandas; print(f'pandas: {pandas.__version__}')" >> "$CONFIG_LOG" 2>&1 || echo "pandas: not available" >> "$CONFIG_LOG"
python3 -c "import openpyxl; print(f'openpyxl: {openpyxl.__version__}')" >> "$CONFIG_LOG" 2>&1 || echo "openpyxl: not available" >> "$CONFIG_LOG"
python3 -c "import gtfparse; print(f'gtfparse: {gtfparse.__version__}')" >> "$CONFIG_LOG" 2>&1 || echo "gtfparse: not available" >> "$CONFIG_LOG"
echo "" >> "$CONFIG_LOG"

# Step 1: Find and copy input BED files
echo "Step 1: Finding and copying input BED files..."
BED_FILES=($(find "$BED_DIR" -name "*.bed" -type f | grep -v "/current_used_file/" | sort))

if [ ${#BED_FILES[@]} -eq 0 ]; then
    echo "Error: No BED files found in $BED_DIR"
    exit 1
fi

# Copy input BED files to output directory
echo "Copying ${#BED_FILES[@]} BED files to input_files directory..."
for file in "${BED_FILES[@]}"; do
    cp "$file" "$INPUT_COPY_DIR/"
    echo "  - Copied $(basename "$file")"
done

# Copy GTF file as well
cp "$GTF_FILE" "$INPUT_COPY_DIR/"
echo "  - Copied $(basename "$GTF_FILE")"

# Copy all necessary scripts for reproducibility
echo "Copying scripts for reproducibility..."
SCRIPT_DIR="$(dirname "$0")"
cp "$SCRIPT_DIR/bed_analysis_pipeline.sh" "$SCRIPTS_COPY_DIR/"
cp "$SCRIPT_DIR/merge_bed_files.sh" "$SCRIPTS_COPY_DIR/"
cp "$SCRIPT_DIR/generate_bed_excel_report.py" "$SCRIPTS_COPY_DIR/"
echo "  - Copied bed_analysis_pipeline.sh"
echo "  - Copied merge_bed_files.sh"
echo "  - Copied generate_bed_excel_report.py"

# Log input files information
echo "=== Input Files ===" >> "$CONFIG_LOG"
echo "Total BED files found: ${#BED_FILES[@]}" >> "$CONFIG_LOG"
echo "BED files:" >> "$CONFIG_LOG"
for file in "${BED_FILES[@]}"; do
    file_size=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null || echo "unknown")
    line_count=$(wc -l < "$file")
    echo "  - $(basename "$file") (${line_count} lines, ${file_size} bytes)" >> "$CONFIG_LOG"
done
echo "GTF file: $(basename "$GTF_FILE")" >> "$CONFIG_LOG"
gtf_size=$(stat -f%z "$GTF_FILE" 2>/dev/null || stat -c%s "$GTF_FILE" 2>/dev/null || echo "unknown")
gtf_lines=$(wc -l < "$GTF_FILE")
echo "  - Size: ${gtf_lines} lines, ${gtf_size} bytes" >> "$CONFIG_LOG"
echo "" >> "$CONFIG_LOG"

# Step 2: Merge BED files
echo "Step 2: Merging BED files..."

echo "Found ${#BED_FILES[@]} BED files:"
for file in "${BED_FILES[@]}"; do
    echo "  - $(basename "$file")"
done
echo ""

# Run merge script
SCRIPT_DIR="$(dirname "$0")"
if [ ! -f "$SCRIPT_DIR/merge_bed_files.sh" ]; then
    echo "Error: merge_bed_files.sh not found in $SCRIPT_DIR"
    exit 1
fi

echo "Running merge script..."
bash "$SCRIPT_DIR/merge_bed_files.sh" "${BED_FILES[@]}" "${MERGED_BED%.*}"

# Log merge results
echo "=== Merge Results ===" >> "$CONFIG_LOG"
if [ -f "$MERGED_LOG" ]; then
    cat "$MERGED_LOG" >> "$CONFIG_LOG"
else
    echo "Merge log not found" >> "$CONFIG_LOG"
fi
echo "" >> "$CONFIG_LOG"
echo ""

# Step 3: Generate comprehensive Excel report
echo "Step 3: Generating comprehensive Excel report with statistics..."

# Check if required Python packages are available
python3 -c "import pandas, openpyxl" 2>/dev/null || {
    echo "Warning: pandas and/or openpyxl not available. Installing..."
    pip3 install pandas openpyxl --user
}

# Check for gtfparse (optional but recommended)
python3 -c "import gtfparse" 2>/dev/null || {
    echo "Note: gtfparse not available. Will use basic GTF parsing."
    echo "For better GTF handling, install with: pip install gtfparse"
}

# Prepare BED files list (include merged file)
ALL_BED_FILES=("${BED_FILES[@]}" "$MERGED_BED")

# Generate Excel file using the separate Python script
SCRIPT_DIR="$(dirname "$0")"
python3 "$SCRIPT_DIR/generate_bed_excel_report.py" \
    "${ALL_BED_FILES[@]}" \
    --gtf "$GTF_FILE" \
    --output "$EXCEL_OUTPUT" \
    --coding-bed "$CODING_BED"

# Log final results
echo "=== Final Results ===" >> "$CONFIG_LOG"
echo "Analysis completed on: $(date)" >> "$CONFIG_LOG"
echo "Processing time: $(($(date +%s) - start_time)) seconds" >> "$CONFIG_LOG"
echo "Total BED files processed: ${#BED_FILES[@]}" >> "$CONFIG_LOG"
if [ -f "$MERGED_BED" ]; then
    echo "Merged BED file peaks: $(wc -l < "$MERGED_BED")" >> "$CONFIG_LOG"
fi
echo "" >> "$CONFIG_LOG"

echo "=== Output Directory Structure ===" >> "$CONFIG_LOG"
echo "Results directory: $RESULTS_DIR" >> "$CONFIG_LOG"
echo "Input files copy: $INPUT_COPY_DIR" >> "$CONFIG_LOG"
echo "Scripts copy: $SCRIPTS_COPY_DIR" >> "$CONFIG_LOG"
echo "Configuration log: $CONFIG_LOG" >> "$CONFIG_LOG"
echo "" >> "$CONFIG_LOG"

echo "=== Reproduction Instructions ===" >> "$CONFIG_LOG"
echo "To reproduce this analysis, run the following command from the output directory:" >> "$CONFIG_LOG"
echo "" >> "$CONFIG_LOG"
echo "chmod +x scripts/bed_analysis_pipeline.sh" >> "$CONFIG_LOG"
echo "./scripts/bed_analysis_pipeline.sh input_files input_files/$(basename "$GTF_FILE") reproduced_output" >> "$CONFIG_LOG"
echo "" >> "$CONFIG_LOG"

# Clean up temporary files
rm -rf "$TEMP_DIR"

echo ""
echo "=== Pipeline Completed Successfully ==="
echo ""
echo "📁 Output Directory Structure:"
echo "   $OUTPUT_DIR/"
echo "   ├── analysis_config.log"
echo "   ├── results/"
echo "   │   ├── coding_regions.bed"
echo "   │   ├── merged_bed_files.bed"
echo "   │   ├── merged_bed_files.log"
echo "   │   └── bed_analysis_report.xlsx"
echo "   ├── input_files/"
echo "   │   ├── [copied BED files]"
echo "   │   └── [copied GTF file]"
echo "   └── scripts/"
echo "       ├── bed_analysis_pipeline.sh"
echo "       ├── merge_bed_files.sh"
echo "       └── generate_bed_excel_report.py"
echo ""

# Display summary statistics
echo "=== Summary Statistics ==="
echo "Total BED files processed: ${#BED_FILES[@]}"
if [ -f "$MERGED_BED" ]; then
    echo "Merged BED file peaks: $(wc -l < "$MERGED_BED")"
fi
echo ""
echo "📊 Detailed analysis available in: $EXCEL_OUTPUT"
echo "📋 Complete configuration log: $CONFIG_LOG"
echo "📂 Input files backed up in: $INPUT_COPY_DIR"
echo "🔧 Scripts backed up in: $SCRIPTS_COPY_DIR"
echo "🔄 Reproduction command in log file (final section)"
