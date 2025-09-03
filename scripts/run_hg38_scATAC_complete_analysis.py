#!/usr/bin/env python3
"""
Complete Human hg38 scATAC Analysis Pipeline

Generates the same comprehensive output as the main bed_analysis_pipeline.sh but
specifically for human hg38 scATAC data with existing merged file.

Creates a complete reproducible analysis package with:
- Individual Excel sheets for each BED file
- Statistics summary with coding region overlap analysis
- All input files copied to output directory
- All scripts copied for reproducibility
- Configuration log with reproduction instructions
"""

import os
import sys
import shutil
import time
import subprocess
from pathlib import Path
from datetime import datetime


def ensure_path_in_syspath(path: Path) -> None:
    """Add path to sys.path if not already present"""
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))


def get_software_versions():
    """Get versions of all software used"""
    versions = {}
    
    # Python version
    versions['python'] = f"Python {sys.version.split()[0]}"
    
    # Try to get bedtools version
    try:
        result = subprocess.run(['bedtools', '--version'], capture_output=True, text=True)
        versions['bedtools'] = result.stdout.strip() if result.returncode == 0 else "not available"
    except:
        versions['bedtools'] = "not available"
    
    # Try to get package versions
    try:
        import pandas
        versions['pandas'] = f"pandas: {pandas.__version__}"
    except:
        versions['pandas'] = "pandas: not available"
    
    try:
        import openpyxl
        versions['openpyxl'] = f"openpyxl: {openpyxl.__version__}"
    except:
        versions['openpyxl'] = "openpyxl: not available"
    
    try:
        import gtfparse
        versions['gtfparse'] = f"gtfparse: {gtfparse.__version__}"
    except:
        versions['gtfparse'] = "gtfparse: not available"
    
    return versions


def copy_input_files(input_files, input_copy_dir):
    """Copy all input files to the input directory"""
    print("Copying input files...")
    input_copy_dir.mkdir(parents=True, exist_ok=True)
    
    copied_files = []
    for file_path in input_files:
        if file_path.exists():
            dest = input_copy_dir / file_path.name
            shutil.copy2(file_path, dest)
            copied_files.append(file_path.name)
            print(f"  - Copied {file_path.name}")
        else:
            print(f"  - Warning: {file_path} not found, skipping")
    
    return copied_files


def copy_scripts(scripts_copy_dir):
    """Copy all necessary scripts for reproducibility"""
    print("Copying scripts for reproducibility...")
    scripts_copy_dir.mkdir(parents=True, exist_ok=True)
    
    script_dir = Path(__file__).resolve().parent
    scripts_to_copy = [
        "generate_bed_excel_report.py",
        "run_hg38_scATAC_complete_analysis.py"
    ]
    
    copied_scripts = []
    for script_name in scripts_to_copy:
        script_path = script_dir / script_name
        if script_path.exists():
            dest = scripts_copy_dir / script_name
            shutil.copy2(script_path, dest)
            copied_scripts.append(script_name)
            print(f"  - Copied {script_name}")
    
    return copied_scripts


def get_file_stats(file_path):
    """Get file statistics (size, line count)"""
    try:
        stat = file_path.stat()
        size = stat.st_size
        
        with open(file_path, 'r') as f:
            line_count = sum(1 for _ in f)
        
        return size, line_count
    except:
        return 0, 0


def create_config_log(config_log_path, start_time, bed_files, gtf_file, merged_bed, 
                     output_dir, results_dir, input_copy_dir, scripts_copy_dir):
    """Create comprehensive configuration log"""
    
    versions = get_software_versions()
    
    with open(config_log_path, 'w') as f:
        f.write("=== Human hg38 scATAC Analysis Configuration Log ===\n")
        f.write(f"Generated on: {datetime.now().strftime('%a %b %d %I:%M:%S %p %Z %Y')}\n")
        f.write("Pipeline version: 1.0 (Human hg38 scATAC specific)\n")
        f.write(f"User: {os.getenv('USER', 'unknown')}\n")
        f.write(f"Host: {os.uname().nodename}\n")
        f.write(f"Working directory: {os.getcwd()}\n")
        f.write("\n")
        
        f.write("=== Input Configuration ===\n")
        f.write(f"BED files directory: bed_file_merger/human\n")
        f.write(f"GTF file: {gtf_file}\n")
        f.write(f"Output directory: {output_dir}\n")
        f.write(f"Merged BED file: {merged_bed}\n")
        f.write("\n")
        
        f.write("=== Software Versions ===\n")
        for name, version in versions.items():
            f.write(f"{version}\n")
        f.write("\n")
        
        f.write("=== Input Files ===\n")
        f.write(f"Total BED files: {len(bed_files)}\n")
        f.write("BED files:\n")
        for bed_file in bed_files:
            if bed_file.exists():
                size, lines = get_file_stats(bed_file)
                f.write(f"  - {bed_file.name} ({lines} lines, {size} bytes)\n")
            else:
                f.write(f"  - {bed_file.name} (not found)\n")
        
        if gtf_file.exists():
            gtf_size, gtf_lines = get_file_stats(gtf_file)
            f.write(f"GTF file: {gtf_file.name}\n")
            f.write(f"  - Size: {gtf_lines} lines, {gtf_size} bytes\n")
        f.write("\n")
        
        f.write("=== Analysis Results ===\n")
        f.write(f"Analysis completed on: {datetime.now().strftime('%a %b %d %I:%M:%S %p %Z %Y')}\n")
        f.write(f"Processing time: {int(time.time() - start_time)} seconds\n")
        f.write(f"Total BED files processed: {len(bed_files)}\n")
        
        if merged_bed.exists():
            _, merged_lines = get_file_stats(merged_bed)
            f.write(f"Merged BED file peaks: {merged_lines}\n")
        f.write("\n")
        
        f.write("=== Output Directory Structure ===\n")
        f.write(f"Results directory: {results_dir}\n")
        f.write(f"Input files copy: {input_copy_dir}\n")
        f.write(f"Scripts copy: {scripts_copy_dir}\n")
        f.write(f"Configuration log: {config_log_path}\n")
        f.write("\n")
        
        f.write("=== Reproduction Instructions ===\n")
        f.write("To reproduce this analysis, run the following command from the output directory:\n")
        f.write("\n")
        f.write("# Activate virtual environment if needed\n")
        f.write("# source /path/to/your/venv/bin/activate\n")
        f.write("\n")
        f.write("chmod +x scripts/run_hg38_scATAC_complete_analysis.py\n")
        f.write("python3 scripts/run_hg38_scATAC_complete_analysis.py --input-dir input_files --output-dir reproduced_output\n")
        f.write("\n")


def main() -> int:
    import argparse
    
    parser = argparse.ArgumentParser(description='Complete Human hg38 scATAC Analysis Pipeline')
    parser.add_argument('--input-dir', default='bed_file_merger/human',
                       help='Directory containing input BED and GTF files')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for complete analysis package')
    
    args = parser.parse_args()
    
    start_time = time.time()
    
    # Setup paths
    if Path(args.input_dir).is_absolute():
        input_dir = Path(args.input_dir)
    else:
        project_root = Path("/home/ls/toozig/gonen-lab/users/toozig/projects/dsd-viz")
        input_dir = project_root / args.input_dir
    
    output_dir = Path(args.output_dir).resolve()
    
    # Define input files
    merged_bed = input_dir / "human_hg38_scATAC_all_combined_fix_id.bed"
    celltype_beds = [
        input_dir / "human_hg38_scATAC_seq_peaks_Sertoli.bed",
        input_dir / "human_hg38_scATAC_seq_peaks_XX_early_supporting.bed",
        input_dir / "human_hg38_scATAC_seq_peaks_XY_early_supporting.bed",
        input_dir / "human_hg38_scATAC_seq_peaks_pre_granulosa.bed",
    ]
    gtf_file = input_dir / "ncbiRefSeq_hg38_23_11_2024.gtf"
    
    # All BED files to analyze (merged + individual cell types)
    all_bed_files = [merged_bed] + celltype_beds
    all_input_files = all_bed_files + [gtf_file]
    
    # Validate inputs
    missing_files = [f for f in all_input_files if not f.exists()]
    if missing_files:
        print("Error: Missing required input files:")
        for f in missing_files:
            print(f"  - {f}")
        return 1
    
    # Create output directory structure
    results_dir = output_dir / "results"
    input_copy_dir = output_dir / "input_files"
    scripts_copy_dir = output_dir / "scripts"
    
    output_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    
    print("=== Human hg38 scATAC Complete Analysis Pipeline ===")
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print("")
    
    # Step 1: Copy input files
    print("Step 1: Copying input files...")
    copied_files = copy_input_files(all_input_files, input_copy_dir)
    
    # Step 2: Copy scripts for reproducibility
    copied_scripts = copy_scripts(scripts_copy_dir)
    
    # Step 3: Generate Excel report
    print("Step 3: Generating comprehensive Excel report with coding statistics...")
    
    # Import the Excel generation functions
    scripts_dir = Path(__file__).resolve().parent
    ensure_path_in_syspath(scripts_dir)
    
    try:
        import generate_bed_excel_report as report
    except Exception as e:
        print(f"Error: failed to import generate_bed_excel_report: {e}")
        return 1
    
    # Setup output files
    excel_output = results_dir / "hg38_scATAC_analysis_report.xlsx"
    coding_bed = results_dir / "hg38_coding_regions_refseq.bed"
    
    # Extract coding regions if not present
    if not coding_bed.exists():
        try:
            report.extract_coding_regions_from_gtf(str(gtf_file), str(coding_bed))
        except Exception as e:
            print(f"Error while extracting coding regions: {e}")
            return 1
    
    # Generate Excel report
    try:
        bed_files_str = [str(f) for f in all_bed_files]
        report.generate_excel_report(bed_files_str, str(coding_bed), str(excel_output))
    except Exception as e:
        print(f"Error while generating Excel report: {e}")
        return 1
    
    # Step 4: Create configuration log
    config_log = output_dir / "analysis_config.log"
    create_config_log(config_log, start_time, all_bed_files, gtf_file, merged_bed,
                     output_dir, results_dir, input_copy_dir, scripts_copy_dir)
    
    print("")
    print("=== Pipeline Completed Successfully ===")
    print("")
    print("📁 Output Directory Structure:")
    print(f"   {output_dir}/")
    print("   ├── analysis_config.log")
    print("   ├── results/")
    print("   │   ├── hg38_coding_regions_refseq.bed")
    print("   │   └── hg38_scATAC_analysis_report.xlsx")
    print("   ├── input_files/")
    print("   │   ├── [copied BED files]")
    print("   │   └── [copied GTF file]")
    print("   └── scripts/")
    print("       ├── generate_bed_excel_report.py")
    print("       └── run_hg38_scATAC_complete_analysis.py")
    print("")
    
    print("=== Summary Statistics ===")
    print(f"Total BED files processed: {len(all_bed_files)}")
    if merged_bed.exists():
        _, merged_lines = get_file_stats(merged_bed)
        print(f"Merged BED file peaks: {merged_lines}")
    print("")
    print(f"📊 Detailed analysis available in: {excel_output}")
    print(f"📋 Complete configuration log: {config_log}")
    print(f"📂 Input files backed up in: {input_copy_dir}")
    print(f"🔧 Scripts backed up in: {scripts_copy_dir}")
    print(f"🔄 Reproduction command in log file (final section)")
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


