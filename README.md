# BED File Merger & Analysis Tool

A comprehensive Python package for merging, analyzing, and generating reports from BED files with coding region analysis. This tool was developed for genomics research and provides standardized analysis workflows for ATAC-seq and other genomic interval data.

## 🚀 Installation & Usage

### Quick Start

```bash
# Install in development mode
pip install -e .

# Run analysis with YAML configuration
python -m bed_file_merger.cli run-config --config configs/your_config.yaml
```

### Requirements

- Python 3.9+
- bedtools (for merging operations)

**Dependencies** (automatically installed):
- pandas ≥2.0.0, numpy <2, openpyxl ≥3.1.0
- pydantic ≥2.5.0, typer ≥0.12.0  
- gtfparse ≥2.5.0, pybedtools ≥0.9.1

## 📋 Configuration

Define your analysis using YAML configuration files:

```yaml
input_config:
  bed_files:
    - path: /path/to/sample1.bed
      name: "Sample 1"  
      extra_columns: ['ID']         # Optional: columns beyond chr, start, end
    - path: /path/to/sample2.bed
      name: "Sample 2"
      extra_columns: ['peak_score', 'ID']
  refseq_gtf: /path/to/refseq.gtf   # GTF file for coding region analysis
  copy_input: true                  # Copy inputs for provenance
  add_id: true                      # Generate unique peak IDs
  id_column_name: id

output_dir: /path/to/output
genome_build: hg38                  # hg38 or mm10
tmp_dir: /tmp

# Optional: Merge overlapping regions using bedtools
merge:
  enabled: true
  out_path: /path/to/output/merged_inputs.bed
  bedtools_opts: ""
```

## 📊 Output

The package generates comprehensive analysis reports:

### 1. **Excel Report** (`report.xlsx`)
- **Individual sheets** for each input BED file with statistics
- **Summary statistics** with coding region overlap analysis  
- **Merged results** sheet (if merging enabled)

### 2. **Merged BED File** (if enabled)
- Combined regions from all input files
- Overlapping intervals merged using bedtools
- Sorted and validated output format

### 3. **Provenance Files** (if `copy_input: true`)
- Input files preserved in `<output_dir>/input_files/`
- Configuration file saved for reproducibility

## 🔧 Available Commands

```bash
# Main analysis command (recommended)
python -m bed_file_merger.cli run-config --config configs/example.yaml

# Extract coding regions from GTF file  
python -m bed_file_merger.cli extract-coding --gtf refseq.gtf --output coding_regions.bed

# Help and command reference
python -m bed_file_merger.cli --help
```

## 📚 Example Workflows

### Human hg38 Analysis
```bash
# Create configuration
cat > configs/human_analysis.yaml << EOF
input_config:
  bed_files:
    - path: data/human_peaks.bed
      name: "Human scATAC"
      extra_columns: ['ID']
  refseq_gtf: data/hg38_refseq.gtf
  copy_input: true
  add_id: true

output_dir: results/human_analysis
genome_build: hg38
merge:
  enabled: true
EOF

# Run analysis
python -m bed_file_merger.cli run-config --config configs/human_analysis.yaml
```

### Cross-Species Analysis
```bash
# Example: Mouse data converted to human coordinates + native human data
cat > configs/cross_species.yaml << EOF
input_config:
  bed_files:
    - path: data/human_native.bed
      name: "Human Native"
      extra_columns: ['ID']
    - path: data/mouse_converted_to_hg38.bed  # Pre-converted using utils
      name: "Mouse Converted"
      extra_columns: ['ID']
  refseq_gtf: data/hg38_refseq.gtf

output_dir: results/cross_species_analysis  
genome_build: hg38
merge:
  enabled: true
EOF

python -m bed_file_merger.cli run-config --config configs/cross_species.yaml
```

## 🛠️ Utility Scripts

This package includes specialized utility scripts for preprocessing and analysis tasks. See **[`utils/README.md`](utils/README.md)** for detailed documentation on:

- **`mm10_to_hg38_liftover.sh`** - Cross-species coordinate conversion
- **`normalize_mm10_beds.sh`** - BED file preprocessing and normalization  
- **`compare_bed_files.sh`** - Detailed comparison between BED files

These utilities complement the main package for complete genomics workflows.

## 🏗️ Project Structure

```
bed_file_merger/            # Main Python package
├── cli.py                  # Command-line interface
├── models.py              # Data models and validation
├── stats.py               # Statistical analysis functions
├── report.py              # Excel report generation
├── coding.py              # Coding region analysis
└── ...

configs/                    # YAML configuration examples
├── hg38_human.yaml
├── mm10_mice.yaml  
└── merging_example.yaml

utils/                      # Utility scripts (see utils/README.md)
├── mm10_to_hg38_liftover.sh
├── normalize_mm10_beds.sh
└── compare_bed_files.sh

results/                    # Analysis outputs
data/                       # Input data files
```

## 🔬 Research Applications

This package was developed for genomics research and is particularly suited for:

- **ATAC-seq peak analysis** with coding region overlap assessment
- **Cross-species genomic comparisons** (mouse ↔ human)
- **Multi-sample peak merging** with provenance tracking
- **Standardized reporting** for genomics publications

## ⚠️ Important Notes

- **Use YAML configurations** - this is the supported and recommended approach
- **Absolute paths recommended** in YAML configs to avoid resolution issues
- **bedtools required** for merge operations (`conda install bedtools`)
- **RefSeq GTF format** expected for coding region analysis
- All coordinate sorting and validation handled automatically

## 🐛 Troubleshooting

**Python module not found:**
```bash
# Install in development mode
pip install -e .
```

**bedtools not found:**
```bash
# Install via conda
conda install bedtools
# Or via system package manager
sudo apt install bedtools  # Ubuntu/Debian
```

**Empty output or errors:**
- Verify input BED files exist and contain data
- Check file paths in YAML configuration
- Ensure input files are tab-separated with chr/start/end columns
- Use absolute paths in YAML configurations

---

**For utility script documentation, see [`utils/README.md`](utils/README.md)**

*This package was developed for genomics research applications and follows best practices for reproducible computational biology.*