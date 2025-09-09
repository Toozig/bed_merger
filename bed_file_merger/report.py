from __future__ import annotations

"""Excel report generation for BED analysis."""

from pathlib import Path
from typing import Iterable, List, Optional

import pandas as pd

ROUND_DECIMALS = 2

# Column names (centralized for readability)
COL_DISPLAY_NAME = "Display Name"
COL_FILE_NAME = "File Name"
COL_TOTAL_ACCESSIBLE_BP = "Total accessible BP"
COL_MIN_PEAK_LENGTH = "Minimal peak Length"
COL_MAX_PEAK_LENGTH = "Maximal peak Length"
COL_MEAN_PEAK_LENGTH = "Mean peak Length"
COL_MEDIAN_PEAK_LENGTH = "Median peak Length"
COL_ACCESSIBLE_CODING_PEAKS = "# Accessible coding Peaks"
COL_ACCESSIBLE_CODING_BP = "Accessible coding BP"
COL_ACCESSIBLE_CODING_PCT = "% Accessible coding\n( Accessible coding BP / Total accessible BP )"
COL_ACCESSIBLE_GENOME_PCT = "% Accessible genome\n( Total accessible BP / Total genome BP )"
COL_CODING_OF_ACCESSIBLE_GENOME_PCT = "% Coding accessible genome (BP)\n(Accessible coding BP / Total genome BP)"
COL_NON_ACCESSIBLE_CODING_PCT_TMPL = "% Non-accessible coding (out of total  {build} bp)\n(Non-accessible coding BP / Total genome BP)"
COL_TOTAL_CODING_PCT_TMPL = "% Total coding BP (out of total {build} bp)"
COL_TOTAL_GENOME_BP_TMPL = "Total {build} bp"
SHEET_SUMMARY_NAME = "Statistics_Summary"
SHEET_DICTIONARY_NAME = "Column_Dictionary"

# Column descriptions (human-readable)
COL_DESCRIPTIONS_STATIC = {
    COL_DISPLAY_NAME: "Configured display name for each input BED.",
    COL_FILE_NAME: "Basename of the original BED file.",
    "Peaks": "Number of peaks (rows) in the BED file.",
    COL_TOTAL_ACCESSIBLE_BP: "Sum of peak lengths (bp) across all peaks.",
    COL_MIN_PEAK_LENGTH: "Minimum peak length in bp.",
    COL_MAX_PEAK_LENGTH: "Maximum peak length in bp.",
    COL_MEAN_PEAK_LENGTH: "Mean peak length in bp.",
    COL_MEDIAN_PEAK_LENGTH: "Median peak length in bp.",
    COL_ACCESSIBLE_CODING_PEAKS: "Number of peaks that overlap coding regions (any overlap).",
    COL_ACCESSIBLE_CODING_BP: "Total base pairs within peaks that overlap coding regions (clipped to the overlap).",
    COL_ACCESSIBLE_CODING_PCT: "Accessible coding BP divided by Total accessible BP (×100).",
    COL_ACCESSIBLE_GENOME_PCT: "Total accessible BP divided by Total genome BP (×100).",
    COL_CODING_OF_ACCESSIBLE_GENOME_PCT: "Accessible coding BP divided by Total genome BP (×100).",
}

def _format_genome_dependent_descriptions(genome_build: str) -> dict:
    return {
        COL_NON_ACCESSIBLE_CODING_PCT_TMPL.format(build=genome_build): (
            "Non-accessible coding BP (Total coding BP − Accessible coding BP) "
            "divided by Total genome BP (×100)."
        ),
        COL_TOTAL_CODING_PCT_TMPL.format(build=genome_build): (
            "Total coding BP divided by Total genome BP (×100)."
        ),
        COL_TOTAL_GENOME_BP_TMPL.format(build=genome_build): (
            "Total genome size in base pairs for the specified genome build."
        ),
    }

def _build_column_dictionary(summary_df: pd.DataFrame, genome_build: str) -> pd.DataFrame:
    """Create a 2-column DataFrame mapping column names to descriptions.

    Only includes the columns present in the summary DataFrame, preserving order.
    """
    mapping = {**COL_DESCRIPTIONS_STATIC, **_format_genome_dependent_descriptions(genome_build)}
    rows = []
    for col in list(summary_df.columns):
        desc = mapping.get(col, "")
        rows.append({"Column": col, "Description": desc})
    return pd.DataFrame(rows, columns=["Column", "Description"])

def process_per_file_frames(per_file_dfs: Iterable[pd.DataFrame]) -> List[pd.DataFrame]:
    """Return processed per-file DataFrames (add size column only)."""
    processed: List[pd.DataFrame] = []
    for df in per_file_dfs:
        tmp = df.copy()
        if "chr" not in tmp.columns and not tmp.empty:
            tmp.columns = [
                "chr", "start", "end",
                *[f"col_{i}" for i in range(4, len(tmp.columns) + 1)],
            ][: len(tmp.columns)]
        if "start" in tmp.columns and "end" in tmp.columns:
            tmp.insert(3, "size", (tmp["end"] - tmp["start"]).astype(int))
        processed.append(tmp)
    return processed


def compute_summary_df(stats_df: pd.DataFrame,
                       genome_total_bp: int,
                       coding_total_bp: int,
                       genome_build: str,
                       display_names: Optional[List[str]]
                       ) -> pd.DataFrame:

    preferred = [
        "file_name",
        "n_peaks",
        "total_bp",
        "min_length",
        "max_length",
        "mean_length",
        "median_length",
        "coding_peaks",
        "coding_bp",
        "coding_percent",
    ]
    cols = [c for c in preferred if c in stats_df.columns]
    summary = stats_df[cols].copy()
    if display_names is not None:
        summary.insert(0, COL_DISPLAY_NAME, display_names[: len(summary)])

    summary.rename(columns={
        "file_name": COL_FILE_NAME,
        "n_peaks": "Peaks",
        "total_bp": COL_TOTAL_ACCESSIBLE_BP,
        "min_length": COL_MIN_PEAK_LENGTH,
        "max_length": COL_MAX_PEAK_LENGTH,
        "mean_length": COL_MEAN_PEAK_LENGTH,
        "median_length": COL_MEDIAN_PEAK_LENGTH,
        "coding_peaks": COL_ACCESSIBLE_CODING_PEAKS,
        "coding_bp": COL_ACCESSIBLE_CODING_BP,
        "coding_percent": COL_ACCESSIBLE_CODING_PCT,
    }, inplace=True)

    # Genome-derived metrics
    
    summary[COL_ACCESSIBLE_GENOME_PCT] = (
        summary[COL_TOTAL_ACCESSIBLE_BP] / genome_total_bp * 100
    ).round(ROUND_DECIMALS)

    summary[COL_CODING_OF_ACCESSIBLE_GENOME_PCT] = (
        summary[COL_ACCESSIBLE_CODING_BP] / genome_total_bp * 100
    ).round(ROUND_DECIMALS)

    summary[COL_NON_ACCESSIBLE_CODING_PCT_TMPL.format(build=genome_build)] = (
        (coding_total_bp - summary[COL_ACCESSIBLE_CODING_BP]) / genome_total_bp * 100
    ).round(ROUND_DECIMALS)
    summary[COL_TOTAL_CODING_PCT_TMPL.format(build=genome_build)] = round(
        coding_total_bp / genome_total_bp * 100, ROUND_DECIMALS
    )
    summary[COL_TOTAL_GENOME_BP_TMPL.format(build=genome_build)] = genome_total_bp

    # print the summary as markdown table
    print(summary.to_string())

    return summary


def save_excel_report(processed_dfs: Iterable[pd.DataFrame],
                      sheet_names: Iterable[str],
                      summary_df: pd.DataFrame,
                      output_excel: Path) -> None:
    """Write per-file sheets and the provided summary to an Excel file.

    Additionally writes a CSV copy of the summary (raw numbers, no formatting)
    next to the Excel file with suffix "_Statistics_Summary.csv".
    """
    output_excel = Path(output_excel)
    output_excel.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
        for df, name in zip(processed_dfs, sheet_names):
            sheet_name = name.replace(".bed", "")[:31]
            df.to_excel(writer, sheet_name=sheet_name, index=False)
        if not summary_df.empty:
            # Write summary to Excel
            summary_df.to_excel(writer, sheet_name=SHEET_SUMMARY_NAME, index=False)

            # Apply thousands separator formatting to numeric columns in summary sheet
            ws = writer.sheets.get(SHEET_SUMMARY_NAME)
            if ws is not None:
                # Determine dtypes from DataFrame to choose number format
                for col_idx, col_name in enumerate(summary_df.columns, start=1):
                    dtype = summary_df[col_name].dtype
                    if pd.api.types.is_integer_dtype(dtype):
                        num_fmt = '#,##0'
                    elif pd.api.types.is_float_dtype(dtype):
                        num_fmt = '#,##0.00'
                    else:
                        num_fmt = None
                    if num_fmt is not None:
                        # Apply format from row 2 (skip header) to last row
                        for row_idx in range(2, ws.max_row + 1):
                            cell = ws.cell(row=row_idx, column=col_idx)
                            cell.number_format = num_fmt

            # Also write a column dictionary tab
            # Try to infer genome_build string by parsing the last columns that include it; fallback to 'genome'
            build = "genome"
            for col in summary_df.columns[::-1]:
                if col.startswith("Total ") and col.endswith(" bp") and "  " not in col:
                    # e.g., "Total hg38 bp" -> extract hg38
                    try:
                        build = col[len("Total ") : -len(" bp")]
                        break
                    except Exception:
                        pass
            dictionary_df = _build_column_dictionary(summary_df, build)
            dictionary_df.to_excel(writer, sheet_name=SHEET_DICTIONARY_NAME, index=False)

    # Write CSV copy of the summary (raw numbers, no Excel formatting)
    if not summary_df.empty:
        csv_path = output_excel.parent / f"{output_excel.stem}_Statistics_Summary.csv"
        summary_df.to_csv(csv_path, index=False)
