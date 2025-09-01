#!/usr/bin/env python3
"""
Run coding vs non-coding statistics XLSX report for human hg38 scATAC data
using the existing generate_bed_excel_report functions.

Inputs are hardcoded to the current project layout:
- Fixed-ID combined BED
- Per-cell-type scATAC BEDs (dated versions)
- RefSeq hg38 GTF provided by the user

Outputs:
- XLSX summary with per-file sheets and Statistics_Summary sheet
- Reusable coding regions BED derived from the GTF
"""

import os
import sys
from pathlib import Path


def ensure_path_in_syspath(path: Path) -> None:
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))


def main() -> int:
    project_root = Path("/home/ls/toozig/gonen-lab/users/toozig/projects/dsd-viz")

    # Inputs
    fixed_id_dir = project_root / "ido_main_data_dir/human_hg38_scATAC_all_combined_fix_id"
    fixed_id_bed = fixed_id_dir / "human_hg38_scATAC_all_combined_fix_id_20_11_2024.bed"

    # Use provided human directory with cleaned BEDs
    human_dir = project_root / "bed_file_merger/human"
    celltype_beds = [
        human_dir / "human_hg38_scATAC_seq_peaks_Sertoli.bed",
        human_dir / "human_hg38_scATAC_seq_peaks_XX_early_supporting.bed",
        human_dir / "human_hg38_scATAC_seq_peaks_XY_early_supporting.bed",
        human_dir / "human_hg38_scATAC_seq_peaks_pre_granulosa.bed",
    ]

    gtf_file = project_root / "ido_main_data_dir/ncbiRefSeq_hg38/ncbiRefSeq_hg38_23_11_2024.gtf"

    # Outputs
    out_dir = fixed_id_dir / "stats_excel"
    out_dir.mkdir(parents=True, exist_ok=True)
    output_xlsx = out_dir / "hg38_scATAC_coding_stats.xlsx"
    coding_bed = out_dir / "hg38_coding_regions_refseq.bed"

    # Validate inputs early
    missing_inputs = [p for p in [fixed_id_bed, gtf_file] + celltype_beds if not p.exists()]
    if missing_inputs:
        print("Error: Missing required input files:")
        for p in missing_inputs:
            print(f"  - {p}")
        return 1

    # Import functions from generate_bed_excel_report.py in the same scripts dir
    scripts_dir = Path(__file__).resolve().parent
    ensure_path_in_syspath(scripts_dir)
    try:
        import generate_bed_excel_report as report
    except Exception as e:
        print(f"Error: failed to import generate_bed_excel_report: {e}")
        return 1

    # Optional: simple bedtools presence check (script relies on bedtools)
    if os.system("command -v bedtools >/dev/null 2>&1") != 0:
        print("Error: bedtools not found in PATH. Please install bedtools and retry.")
        return 1

    # Build coding regions BED once (if not already present)
    if not coding_bed.exists() or coding_bed.stat().st_size == 0:
        try:
            report.extract_coding_regions_from_gtf(str(gtf_file), str(coding_bed))
        except Exception as e:
            print(f"Error while extracting coding regions: {e}")
            return 1

    # Assemble BED list (fixed-id combined + per-cell-type)
    bed_files = [str(fixed_id_bed)] + [str(p) for p in celltype_beds]

    # Generate Excel report
    try:
        report.generate_excel_report(bed_files, str(coding_bed), str(output_xlsx))
    except Exception as e:
        print(f"Error while generating Excel report: {e}")
        return 1

    print("\nDone.")
    print(f"XLSX: {output_xlsx}")
    print(f"Coding BED: {coding_bed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


