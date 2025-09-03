# BED File Merger - Project Organization

This directory contains all files related to BED file merging and analysis for the mm10 ATAC-seq project.

## Directory Structure

### 📁 `scripts/`
Contains all analysis scripts:
- `merge_bed_files.sh` - Script to merge multiple BED files and generate statistics
- `compare_bed_files.sh` - Script to compare two BED files and generate intersection outputs

### 📁 `merged_outputs/`
Contains merged BED file results:
- `merged_mice_mm10_atac.bed` - Final merged BED file from all input files (130,173 peaks)
- `merged_mice_mm10_atac.log` - Statistics log from the merging process

### 📁 `analysis_results/`
Contains comparison analysis outputs:
- `complete_overlaps_with_ids.bed` - Peaks with 100% complete overlap (83,326 peaks)
- `intersections_with_ids.bed` - All intersecting peaks with IDs (120,911 peaks)
- `merged_only_segments.bed` - Novel segments only in merged file (9,442 peaks)
- `all_segments_combined.bed` - Combined file with all segments (130,353 peaks)
- `comparison_statistics.txt` - Detailed comparison statistics
- `novel_segments_analysis.txt` - Analysis of novel segment characteristics

### 📁 `mice_mm10_atac_seq_bed_Files/`
Contains original input BED files:
- Individual ATAC-seq BED files from different developmental stages
- `current_used_file/mouse_mm10_peaks_withID_1406.bed` - Reference file with IDs

## Key Findings

- **Reference file**: 120,911 peaks (98.9M bp coverage)
- **Merged file**: 130,173 peaks (105.3M bp coverage)
- **100% of reference peaks** are represented in the merged dataset
- **9,442 additional novel peaks** discovered (primarily shorter peaks 201-500bp)
- **Novel segments are legitimate regulatory elements** that were likely filtered out in original analysis

## Usage

To recreate the analysis:

```bash
# Merge BED files
./scripts/merge_bed_files.sh mice_mm10_atac_seq_bed_Files/*.bed output_name

# Compare merged file with reference
./scripts/compare_bed_files.sh merged_file.bed reference_file.bed output_directory
```

## Analysis Date
Generated: $(date)





