# BED File Merger - Project Organization

This directory contains all files related to BED file merging and analysis.

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

To run the YAML-driven CLI:

```bash
python -m bed_file_merger.cli run-config --config configs/hg38_human.yaml
```

### YAML configuration (new schema)

```yaml
input_config:
  bed_files:
    - path: /path/a.bed
      name: A
      extra_columns: [score, id]
    - path: /path/b.bed
      name: B
      extra_columns: null
  refseq_gtf: /path/refseq.gtf
  coding_bed: /path/to/output/hg38_coding_regions_refseq.bed
  copy_input: false

output_dir: /path/to/output
genome_build: hg38

merge:
  enabled: false
  out_path: /path/to/output/merged_inputs.bed
  bedtools_opts: ""

# Names and extra column specifications are provided per file above.
```

Notes:
- `copy_input: true` will copy all inputs to `<output_dir>/input_files/` for provenance; analysis uses original paths.
- `coding_bed` is optional; if omitted, it will be generated under `output_dir`.

## Analysis Date
Generated: $(date)





