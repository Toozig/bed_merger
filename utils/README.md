# Utility Scripts

This directory contains specialized utility scripts that support the main BED File Merger package. These are standalone tools for specific genomic data processing tasks.

## 🔄 `mm10_to_hg38_liftover.sh`

**Purpose:** Converts genomic coordinates from mm10 (mouse) genome to hg38 (human) genome using UCSC's liftOver tool.

**Usage:**
```bash
# Default usage (uses preset file paths)
./utils/mm10_to_hg38_liftover.sh

# Custom usage with all parameters
./utils/mm10_to_hg38_liftover.sh input.bed chain_file.gz output_success.bed output_failed.bed 0.1

# Help and version
./utils/mm10_to_hg38_liftover.sh --help
./utils/mm10_to_hg38_liftover.sh --version
```

**Parameters (all optional):**
1. `input_bed` - Input BED file (default: `results/mm10_merged_peaks_030925.bed`)
2. `chain_file` - Chain file (default: `mm10ToHg38.over.chain.gz`)
3. `output_success` - Successful conversions (default: `results/liftover_mm10_to_hg38_peaks.bed`)
4. `output_failed` - Failed conversions (default: `results/liftover_mm10_unsuccessful.bed`)
5. `min_match` - Minimum match ratio (default: `0.1`)

**Features:**
- Cross-species coordinate conversion
- Multiple mappings per region (`-multiple` flag)
- Conversion statistics reporting
- Input file validation
- Automatic output directory creation

**When to use:** When you need to compare mouse genomic regions with human data by converting mouse coordinates to human reference genome.

---

## 🔧 `normalize_mm10_beds.sh`

**Purpose:** Preprocesses and normalizes mm10 BED files by filtering chromosomal regions and standardizing format.

**Usage:**
```bash
./utils/normalize_mm10_beds.sh
```

**Input:** Reads from `data/mm10_atac_seq/raw_files/*.bed`  
**Output:** Writes to `data/mm10_atac_seq/*_normalized.bed`

**What it does:**
- Filters to keep only chromosomal regions (removes contigs/unmapped sequences starting with `chrUn_`, `chrM`, etc.)
- Retains only the first three columns: `chr`, `start`, `end`
- Standardizes format for downstream analysis
- Processes all BED files in the raw_files directory

**When to use:** As a preprocessing step for raw mm10 ATAC-seq BED files before analysis with the main package.

---

## 📊 `compare_bed_files.sh`

**Purpose:** Performs detailed comparison analysis between two BED files, typically a merged file and a reference file.

**Usage:**
```bash
./utils/compare_bed_files.sh merged.bed reference_with_ids.bed output_directory/
```

**Parameters:**
1. `merged_bed` - Merged BED file to analyze
2. `reference_bed_with_ids` - Reference BED file with ID column
3. `output_dir` - Directory for output files

**Generated outputs:**
- `complete_overlaps_with_ids.bed` - Regions with 100% complete overlap
- `intersections_with_ids.bed` - All intersecting regions with IDs
- `merged_only_segments.bed` - Novel segments only in merged file  
- `all_segments_combined.bed` - Combined file with all segments
- `comparison_statistics.txt` - Detailed comparison statistics
- `novel_segments_analysis.txt` - Analysis of novel segment characteristics

**Features:**
- Complete overlap detection
- Novel segment identification
- Statistical analysis of differences
- ID preservation from reference file
- Comprehensive reporting

**When to use:** When you need to understand the differences between a merged dataset and a reference, particularly to identify novel regions discovered through merging.

---

## 📈 `mm10_to_hg38_analysis.py`

Replicates the `mm10_tohg_38_analysis.ipynb` workflow to produce an Excel report (with an extra tab for unsuccessful conversions) and a filtered BED.

Usage:
```bash
python utils/mm10_to_hg38_analysis.py \
  --original-bed results/mm10_to_hg38_files/mm10_merged_peaks_030925.bed \
  --converted-bed results/mm10_to_hg38_files/liftOveroutput/liftover_mm10_to_hg38_peaks.bed \
  --unsuccessful-bed results/mm10_to_hg38_files/liftOveroutput/liftover_mm10_unsuccseful.bed \
  --output-dir results/mm10_to_hg38_files/liftOverFilteration \
  --rel-tol 0.25 --abs-tol 1000 --min-size-floor 100
```

Outputs:
- `mm10_to_hg38_analysis_results.xlsx`: Tabs `Filtered_Peaks`, `All_Peak`, `All_Statistics`, `Unsuccessful`
- `mm10_to_hg38_filtered.bed`: Filtered peaks (chr, start, end, ID)

## 🔗 Integration with Main Package

These utilities complement the main BED File Merger Python package:

1. **Preprocessing workflow:**
   ```bash
   # 1. Normalize raw mm10 files
   ./utils/normalize_mm10_beds.sh
   
   # 2. Convert to human coordinates (if needed)
   ./utils/mm10_to_hg38_liftover.sh
   
   # 3. Run main analysis
   python -m bed_file_merger.cli run-config --config configs/your_config.yaml
   
   # 4. Compare results (if needed)
   ./utils/compare_bed_files.sh merged_output.bed reference.bed comparison_results/
   ```

2. **Cross-species analysis:**
   - Use `mm10_to_hg38_liftover.sh` to convert mouse data
   - Include converted data in YAML config alongside human data
   - Run integrated analysis with the main package

## Requirements

- **bash** (for shell scripts)
- **bedtools** (for BED operations)
- **liftOver** (UCSC tool, for coordinate conversion)
- Standard Unix tools: `awk`, `grep`, `sort`, `wc`

## Notes

- All scripts include input validation and error handling
- Use absolute paths or ensure you're in the correct working directory
- Scripts are designed to work with the directory structure created by the main package
- For questions about the main analysis package, see the root README.md
