#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

import pandas as pd

# input file
original_bed = Path("results/Supp_data_2_mm10_gnoad_cells_ATAC_seq/Supp_data_2_mm10_gnoad_cells_ATAC_seq.bed")
converted_bed = Path("results/Supp_Data_3_mm10_to_hg38_conversion/lifOver_output/liftover_mm10_to_hg38_peaks.bed")
unsuccessful_bed = Path("results/Supp_Data_3_mm10_to_hg38_conversion/lifOver_output/liftover_mm10_unsuccessful.bed")
output_dir = Path("results/Supp_Data_3_mm10_to_hg38_conversion/liftOver_filteration")

# output file
XL_OUTPUT = output_dir / 'Supp_Data_3_mm10_to_hg38_conversion.xlsx'
BED_FILTERED = output_dir / 'Supp_Data_3_mm10_to_hg38_conversion.bed'


# Defaults (mirroring the notebook constants)
REL_TOL_DEFAULT: float = 0.25   # 25%
ABS_TOL_DEFAULT: int = 1000     # bp
MIN_SIZE_FLOOR_DEFAULT: int = 100


def read_bed(path: Path, names: list[str]) -> pd.DataFrame:
    """Read a BED-like file, skipping comment lines starting with '#'."""
    df = pd.read_csv(path, sep="\t", header=None, names=names, comment="#")
    return df


def add_size(df: pd.DataFrame, size_col: str) -> pd.DataFrame:
    out = df.copy()
    out[size_col] = out["end"].astype(int) - out["start"].astype(int)
    return out


def build_comparison(original_bed: Path, converted_bed: Path) -> pd.DataFrame:
    original = read_bed(original_bed, ["chr", "start", "end", "ID"])
    original = add_size(original, "original_size")

    converted = read_bed(converted_bed, ["chr", "start", "end", "ID", "n_multiple"])
    # Keep IDs unchanged for a correct merge against the original (mm10) IDs
    converted = add_size(converted, "converted_size")

    # Keep sizes and n_multiple
    comp = pd.merge(original[["ID", "original_size"]],
                    converted[["chr", "start", "end", "ID", "n_multiple", "converted_size"]],
                    on="ID",
                    how="right")
    # Robust difference: handle missing sizes by treating them as 0
    orig_size = pd.to_numeric(comp["original_size"], errors="coerce").fillna(0).astype(int)
    conv_size = pd.to_numeric(comp["converted_size"], errors="coerce").fillna(0).astype(int)
    comp["diffrence"] = orig_size - conv_size
    return original, converted, comp


def compute_min_size_threshold(original: pd.DataFrame, rel_tol: float, size_floor: int) -> int:
    minimal_size = int(original["original_size"].min())
    minimal_size = max(int(minimal_size - (minimal_size * rel_tol)), int(size_floor))
    return minimal_size


def apply_filters(comp: pd.DataFrame, abs_tol: int, min_size: int) -> pd.DataFrame:
    def _row_ok(row: pd.Series) -> bool:
        original = int(row["original_size"]) if pd.notna(row["original_size"]) else 0
        converted = int(row["converted_size"]) if pd.notna(row["converted_size"]) else 0
        abs_diff = abs(original - converted)
        abs_diff_filter = abs_diff <= abs_tol
        size_filter = converted >= min_size
        return abs_diff_filter and size_filter

    out = comp.copy()
    out["passed_size_filter"] = out.apply(_row_ok, axis=1)
    out["passed_multiple_filter"] = ~out["ID"].isin(
        out.loc[out["n_multiple"].astype(int) > 1, "ID"]
    )
    return out


def compute_stats(original: pd.DataFrame, converted: pd.DataFrame, comp: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    faild_to_convert = original.shape[0] - converted.shape[0]

    peak_filtering_stats = {
        'metric': ['Total original peaks', 'Total converted peaks', 'Failed to convert'],
        'count': [original.shape[0], converted.shape[0], faild_to_convert],
        'percentage': [100.0, converted.shape[0]/original.shape[0]*100, faild_to_convert/original.shape[0]*100]
    }

    size_filter_passed = comp['passed_size_filter'].sum()
    size_filter_failed = len(comp) - size_filter_passed
    size_filter_stats = {
        'metric': ['Passed size filter', 'Failed size filter'],
        'count': [size_filter_passed, size_filter_failed],
        'percentage': [size_filter_passed/len(comp)*100, size_filter_failed/len(comp)*100]
    }

    multiple_filter_passed = comp['passed_multiple_filter'].sum()
    multiple_filter_failed = len(comp) - multiple_filter_passed
    multiple_filter_stats = {
        'metric': ['Passed multiple mapping filter', 'Failed multiple mapping filter'],
        'count': [multiple_filter_passed, multiple_filter_failed],
        'percentage': [multiple_filter_passed/len(comp)*100, multiple_filter_failed/len(comp)*100]
    }

    both_filters_passed = (comp['passed_size_filter'] & comp['passed_multiple_filter']).sum()
    combined_filter_stats = {
        'metric': ['Passed both filters', 'Failed at least one filter', 'Final conversion success rate'],
        'count': [both_filters_passed, len(comp) - both_filters_passed, both_filters_passed],
        'percentage': [both_filters_passed/len(comp)*100, (len(comp) - both_filters_passed)/len(comp)*100, 
                      both_filters_passed/original.shape[0]*100]
    }

    return (
        pd.DataFrame(peak_filtering_stats),
        pd.DataFrame(size_filter_stats),
        pd.DataFrame(multiple_filter_stats),
        pd.DataFrame(combined_filter_stats),
    )


def write_excel(
    output_xlsx: Path,
    filtered_df: pd.DataFrame,
    original_df: pd.DataFrame,
    converted_df: pd.DataFrame,
    all_stats_df: pd.DataFrame,
    unsuccessful_bed: Path | None,
) -> None:
    """Write Excel with required tabs and formatted statistics.

    Tabs (order):
    - 'filtered peaks': peaks passing filters (first tab)
    - 'original peaks': original mm10 peaks
    - 'converted peaks': successfully converted peaks
    - 'Failed to convert': peaks that failed conversion (BED may contain comment lines '#')
    - 'Statistics': combined statistics with thousands separator formatting
    - 'Explanation': textual description of tabs and filters
    """
    output_xlsx.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_xlsx, engine='openpyxl') as writer:
        # Data tabs (filtered first)
        filt_cols = [c for c in ["chr", "start", "end", "ID", "n_multiple", "converted_size"] if c in filtered_df.columns]
        filtered_df[filt_cols].to_excel(writer, sheet_name='filtered peaks', index=False)

        orig_cols = [c for c in ["chr", "start", "end", "ID"] if c in original_df.columns]
        conv_cols = [c for c in ["chr", "start", "end", "ID", "n_multiple"] if c in converted_df.columns]
        original_df[orig_cols].to_excel(writer, sheet_name='original peaks', index=False)
        converted_df[conv_cols].to_excel(writer, sheet_name='converted peaks', index=False)

        # Failed to convert tab (optional)
        if unsuccessful_bed is not None and Path(unsuccessful_bed).exists():
            df_uns = read_bed(Path(unsuccessful_bed), ["chr", "start", "end", "ID"])  # '#' comments ignored
            df_uns.to_excel(writer, sheet_name='Failed to convert', index=False)

        # Statistics tab with number formatting
        all_stats_df.to_excel(writer, sheet_name='Statistics', index=False)
        ws = writer.sheets.get('Statistics')
        if ws is not None and not all_stats_df.empty:
            for col_idx, col_name in enumerate(all_stats_df.columns, start=1):
                dtype = all_stats_df[col_name].dtype
                if pd.api.types.is_integer_dtype(dtype):
                    num_fmt = '#,##0'
                elif pd.api.types.is_float_dtype(dtype):
                    num_fmt = '#,##0.00'
                else:
                    num_fmt = None
                if num_fmt is not None:
                    for row_idx in range(2, ws.max_row + 1):  # skip header
                        ws.cell(row=row_idx, column=col_idx).number_format = num_fmt

        # Explanation tab (textual)
        explanation_rows = [
            {"Section": "filtered peaks", "Description": "Peaks that passed both filters (size and multi-mapping)."},
            {"Section": "original peaks", "Description": "Original mm10 peaks (chr,start,end,ID)."},
            {"Section": "converted peaks", "Description": "Peaks after liftover to hg38 (includes n_multiple)."},
            {"Section": "Failed to convert", "Description": "Peaks that failed conversion; comment lines ('#') are ignored on load."},
            {"Section": "Statistics", "Description": "Summary counts and percentages; numbers formatted with thousands separators."},
            {"Section": "Filter: size", "Description": "Keeps peaks with converted_size ≥ min_size of smallest original peak; min_size computed from original sizes using relative tolerance (0.25 of original size) and floor."},
            {"Section": "Filter: multiple mapping", "Description": "Removes peaks with n_multiple > 1 (non-unique mappings)."},
        ]
        pd.DataFrame(explanation_rows, columns=["Section", "Description"]).to_excel(
            writer, sheet_name='Explanation', index=False
        )


def get_args():
    p = argparse.ArgumentParser(description='mm10→hg38 conversion analysis (scripted from notebook logic)')
    p.add_argument('--original-bed', default=str(original_bed), help='Original merged peaks BED (chr,start,end,ID)')
    p.add_argument('--converted-bed', default=str(converted_bed), help='Converted peaks BED (chr,start,end,ID,n_multiple)')
    p.add_argument('--unsuccessful-bed', default=str(unsuccessful_bed), help='Unsuccessful converted BED to add as a tab')
    p.add_argument('--output-dir', default=str(output_dir), help='Output directory for analysis results')
    p.add_argument('--rel-tol', type=float, default=REL_TOL_DEFAULT, help=f'Relative tolerance (default {REL_TOL_DEFAULT})')
    p.add_argument('--abs-tol', type=int, default=ABS_TOL_DEFAULT, help=f'Absolute tolerance in bp (default {ABS_TOL_DEFAULT})')
    p.add_argument('--min-size-floor', type=int, default=MIN_SIZE_FLOOR_DEFAULT, help=f'Floor for minimal accepted size (default {MIN_SIZE_FLOOR_DEFAULT})')
    args = p.parse_args()
    return args


def main() -> int:
    args = get_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    original, converted, comp = build_comparison(Path(args.original_bed), Path(args.converted_bed))
    min_size = compute_min_size_threshold(original, args.rel_tol, args.min_size_floor)
    comp = apply_filters(comp, abs_tol=args.abs_tol, min_size=min_size)

    peak_df, size_df, multi_df, both_df = compute_stats(original, converted, comp)
    peak_df['category'] = 'Peak Filtering'
    size_df['category'] = 'Size Filter'
    multi_df['category'] = 'Multiple Mapping Filter'
    both_df['category'] = 'Combined Filter'
    all_stats_df = pd.concat([peak_df, size_df, multi_df, both_df], ignore_index=True)

    # Build filtered/comparison with columns aligned to notebook
    filtered_df = comp[comp['passed_size_filter'] & comp['passed_multiple_filter']].copy()
    filtered_df = filtered_df[['chr', 'start', 'end', 'ID', 'n_multiple', 'converted_size']].copy()
    # rename id
    converted['ID'] = converted['ID'].str.replace('mm10', 'hg38')
    filtered_df['ID'] = filtered_df['ID'].str.replace('mm10', 'hg38')
    # Write outputs
    write_excel(
        XL_OUTPUT,
        filtered_df=filtered_df,
        original_df=original,
        converted_df=converted,
        all_stats_df=all_stats_df,
        unsuccessful_bed=Path(args.unsuccessful_bed) if args.unsuccessful_bed else None,
    )
    filtered_df[['chr', 'start', 'end', 'ID']].to_csv(BED_FILTERED, sep='\t', index=False, header=False)

    print(f"Wrote Excel: {XL_OUTPUT}")
    print(f"Wrote filtered BED: {BED_FILTERED}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())


