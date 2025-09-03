from __future__ import annotations

"""Excel report generation for BED analysis."""

from pathlib import Path
from typing import Iterable, List, Optional

import pandas as pd

ROUND_DECIMALS = 2

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
        summary.insert(0, "Display Name", display_names[: len(summary)])

    summary.rename(columns={
        "file_name": "File Name",
        "n_peaks": "Peaks",
        "total_bp": "Total BP",
        "min_length": "Min Length",
        "max_length": "Max Length",
        "mean_length": "Mean Length",
        "median_length": "Median Length",
        "coding_peaks": "Coding Peaks",
        "coding_bp": "Coding BP",
        "coding_percent": "Coding % accessible",
    }, inplace=True)

    # Genome-derived metrics
    
    summary["Accessible % genome"] = (summary["Total BP"] / 
    float(genome_total_bp) * 100).round(ROUND_DECIMALS)

    summary["% coding of accessible genome "] = (summary["Coding BP"]/ summary["Total BP"] * 100).round(ROUND_DECIMALS)

    summary[f"Non-accessible % coding (out of total  {genome_build} bp)"] = (
        ((coding_total_bp - summary["Coding BP"]) /genome_total_bp ) * 100
    ).round(ROUND_DECIMALS)
    summary[f"Total % coding BP (out of total {genome_build} bp)"] = round(coding_total_bp / genome_total_bp * 100, ROUND_DECIMALS)
    summary[f"Total {genome_build} bp"] = genome_total_bp

    # print the summary as markdown table
    print(summary.to_string())

    return summary


def save_excel_report(processed_dfs: Iterable[pd.DataFrame],
                      sheet_names: Iterable[str],
                      summary_df: pd.DataFrame,
                      output_excel: Path) -> None:
    """Write per-file sheets and the provided summary to an Excel file."""
    output_excel = Path(output_excel)
    output_excel.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
        for df, name in zip(processed_dfs, sheet_names):
            sheet_name = name.replace(".bed", "")[:31]
            df.to_excel(writer, sheet_name=sheet_name, index=False)
        if not summary_df.empty:
            summary_df.to_excel(writer, sheet_name="Statistics_Summary", index=False)
