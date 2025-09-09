"""Plot genome accessibility pie charts from summary CSV.

This script reads a CSV summary table with columns including `Display Name`,
`Accessible coding BP`, `% Total coding BP (out of total hg38 bp)`, and
`Total hg38 bp`. For each selected row (by `Display Name`), it computes three
genome partitions:

1) Non-coding base pairs
2) Inaccessible coding base pairs
3) Accessible coding base pairs

It then renders two variants of a pie chart for each row:
- A labeled pie chart (percentages on wedges) saved as PNG and SVG
- A clear pie chart (no wedge labels, legend only) saved as PNG and SVG

Configuration is controlled via module-level globals at the top of the file.

Usage: Adjust the global variables and run this script with Python.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import math
import os
import re
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import matplotlib.pyplot as plt
import pandas as pd
from pydantic import BaseModel, Field


# ==========================
# Globals (edit as needed)
# ==========================

# Defaults are set via CLI; below are fallbacks used only if CLI not provided
INPUT_CSV: str | None = None
OUTPUT_DIR: str | None = None
DISPLAY_NAMES: List[str] | None = None
USE_LOG: bool = True  # default to log-sizing; can be disabled via --no-log

# Colors used for [non-coding, inaccessible coding, accessible coding]
COLORS = {
    "non_coding_inaccessible": "#8DA0CB",  # Set2 blue
    "non_coding_accessible": "#66C2A5",    # Set2 green
    "coding_inaccessible": "#FC8D62",      # Set2 orange
    "coding_accessible": "#E78AC3",        # Set2 pink
}


class RowModel(BaseModel):
    """Validated representation of a single row in the CSV.

    Attributes
    ----------
    display_name: str
        The display name identifying the row.
    total_genome_bp: int
        Total number of base pairs in the genome (e.g., hg38).
    accessible_coding_bp: int
        Number of base pairs that are coding and accessible.
    total_coding_percent: float
        Percentage of total genome that is coding (0-100).
    """

    display_name: str
    total_genome_bp: int = Field(ge=0)
    accessible_coding_bp: int = Field(ge=0)
    total_accessible_bp: int = Field(ge=0)
    total_coding_percent: float = Field(ge=0.0, le=100.0)

    @property
    def total_coding_bp(self) -> float:
        """Compute total coding base pairs from percent of the genome.

        Returns
        -------
        float
            Total coding base pairs.
        """
        return (self.total_coding_percent / 100.0) * float(self.total_genome_bp)

    @property
    def inaccessible_coding_bp(self) -> float:
        """Compute inaccessible coding base pairs.

        Returns
        -------
        float
            Inaccessible coding base pairs (non-negative).
        """
        value = self.total_coding_bp - float(self.accessible_coding_bp)
        return value if value >= 0 else 0.0

    @property
    def non_coding_bp(self) -> float:
        """Compute non-coding base pairs.

        Returns
        -------
        float
            Non-coding base pairs (non-negative).
        """
        value = float(self.total_genome_bp) - self.total_coding_bp
        return value if value >= 0 else 0.0

    @property
    def accessible_genome_percent(self) -> float:
        """Percentage of the genome that is accessible (0-100)."""
        if self.total_genome_bp == 0:
            return 0.0
        return 100.0 * (float(self.total_accessible_bp) / float(self.total_genome_bp))

    @property
    def coding_accessible_percent(self) -> float:
        """Percentage of the genome that is coding and accessible (0-100)."""
        if self.total_genome_bp == 0:
            return 0.0
        return 100.0 * (float(self.accessible_coding_bp) / float(self.total_genome_bp))

    @property
    def coding_inaccessible_percent(self) -> float:
        """Percentage of the genome that is coding and inaccessible (0-100)."""
        value = self.total_coding_percent - self.coding_accessible_percent
        return value if value >= 0 else 0.0

    @property
    def noncoding_total_percent(self) -> float:
        """Percentage of the genome that is non-coding (0-100)."""
        value = 100.0 - self.total_coding_percent
        return value if value >= 0 else 0.0

    @property
    def noncoding_accessible_percent(self) -> float:
        """Percentage of the genome that is non-coding and accessible (0-100)."""
        value = self.accessible_genome_percent - self.coding_accessible_percent
        return value if value >= 0 else 0.0

    @property
    def noncoding_inaccessible_percent(self) -> float:
        """Percentage of the genome that is non-coding and inaccessible (0-100)."""
        value = self.noncoding_total_percent - self.noncoding_accessible_percent
        return value if value >= 0 else 0.0


def sanitize_display_name(name: str) -> str:
    """Sanitize display name for filenames.

    - Trim whitespace
    - Replace spaces with underscores
    - Remove parentheses and any non [A-Za-z0-9_-]

    Parameters
    ----------
    name : str
        Original display name.

    Returns
    -------
    str
        Sanitized string suitable for filenames.
    """
    trimmed = name.strip()
    underscored = re.sub(r"\s+", "_", trimmed)
    # Remove all but alphanumerics, dash, underscore
    cleaned = re.sub(r"[^A-Za-z0-9_\-]+", "", underscored)
    # Collapse multiple underscores/dashes
    cleaned = re.sub(r"[_\-]{2,}", lambda m: m.group(0)[0], cleaned)
    return cleaned.strip("_- ")


def _normalize_colname(name: str) -> str:
    """Normalize a column name for regex matching: strip quotes, trim, collapse spaces."""
    no_quotes = name.strip().strip('"').strip("'")
    collapsed = re.sub(r"\s+", " ", no_quotes)
    return collapsed


def read_summary_csv(input_csv: str) -> pd.DataFrame:
    """Read the summary CSV and drop spurious non-data rows.

    The file can contain explanatory header lines after the header row.
    We keep only rows where key numeric columns are present and valid.

    Parameters
    ----------
    input_csv : str
        Path to the CSV file.

    Returns
    -------
    pandas.DataFrame
        Cleaned dataframe with consistent columns.
    """
    df = pd.read_csv(input_csv, engine="python")

    # Pre-normalize all column names for easier matching (temporary mapping)
    normalized_map = {c: _normalize_colname(c) for c in df.columns}
    inverse_map = {v: k for k, v in normalized_map.items()}

    # Identify species-specific columns and rename them to canonical names
    # Find Total genome bp column (Total hg38 bp or Total mm10 bp)
    total_genome_norm = next(
        (n for n in normalized_map.values() if re.match(r"^Total (hg38|mm10) bp$", n, flags=re.IGNORECASE)),
        None,
    )
    if total_genome_norm is None:
        raise ValueError("Could not find 'Total <genome> bp' column (hg38/mm10).")
    total_genome_col = inverse_map[total_genome_norm]

    # Find % Total coding BP column
    coding_percent_norm = next(
        (
            n
            for n in normalized_map.values()
            if re.match(r"^% Total coding BP \(out of total (hg38|mm10) bp\)$", n, flags=re.IGNORECASE)
        ),
        None,
    )
    if coding_percent_norm is None:
        raise ValueError(
            "Could not find '% Total coding BP (out of total <genome> bp)' column (hg38/mm10)."
        )
    coding_percent_col = inverse_map[coding_percent_norm]

    # Canonicalize these columns for downstream processing
    df = df.rename(
        columns={
            total_genome_col: "Total genome bp",
            coding_percent_col: "% Total coding BP (out of total genome bp)",
        }
    )

    # Coerce key numeric fields to numeric and filter valid rows
    df["Accessible coding BP"] = pd.to_numeric(df.get("Accessible coding BP"), errors="coerce")
    df["% Total coding BP (out of total genome bp)"] = pd.to_numeric(
        df.get("% Total coding BP (out of total genome bp)"), errors="coerce"
    )
    df["Total genome bp"] = pd.to_numeric(df.get("Total genome bp"), errors="coerce")
    df["Total accessible BP"] = pd.to_numeric(df.get("Total accessible BP"), errors="coerce")

    mask_valid = (
        df["Accessible coding BP"].notna()
        & df["% Total coding BP (out of total genome bp)"].notna()
        & df["Total genome bp"].notna()
        & df["Total accessible BP"].notna()
    )
    df = df.loc[mask_valid].copy()

    # Ensure display name is str and drop empties
    df["Display Name"] = df["Display Name"].astype(str)
    df = df[df["Display Name"].str.strip() != ""]

    return df


def dataframe_to_rows(df: pd.DataFrame, display_names: Sequence[str] | None) -> List[RowModel]:
    """Convert a dataframe to a list of validated RowModel entries.

    Parameters
    ----------
    df : pandas.DataFrame
        The cleaned dataframe.
    display_names : Sequence[str] | None
        Names to select. If None or empty, all rows are used.

    Returns
    -------
    list[RowModel]
        Validated row models.
    """
    if display_names:
        df = df[df["Display Name"].isin(display_names)].copy()

    rows: List[RowModel] = []
    for _, r in df.iterrows():
        model = RowModel(
            display_name=str(r["Display Name"]),
            total_genome_bp=int(r["Total genome bp"]),
            accessible_coding_bp=int(r["Accessible coding BP"]),
            total_accessible_bp=int(r["Total accessible BP"]),
            total_coding_percent=float(r["% Total coding BP (out of total genome bp)"]),
        )
        rows.append(model)
    return rows


def compute_segments(row: RowModel) -> Tuple[List[str], List[float]]:
    """Compute four percentage-based segments for the pie chart.

    Segments (sum to ~100):
    - Non-coding inaccessible
    - Non-coding accessible
    - Coding inaccessible
    - Coding accessible

    Parameters
    ----------
    row : RowModel
        Validated row model.
    Returns
    -------
    tuple[list[str], list[float]]
        Labels and raw percent sizes (0-100) in plotting order.
    """
    labels = [
        "Non-coding inaccessible",
        "Non-coding accessible",
        "Coding inaccessible",
        "Coding accessible",
    ]

    sizes = [
        row.noncoding_inaccessible_percent,
        row.noncoding_accessible_percent,
        row.coding_inaccessible_percent,
        row.coding_accessible_percent,
    ]

    # Ensure no negative tiny values due to floating precision
    sizes = [max(0.0, float(v)) for v in sizes]

    # Normalize minor rounding to exactly 100 if needed
    total = sum(sizes)
    if total > 0 and abs(total - 100.0) > 1e-6:
        sizes = [v * 100.0 / total for v in sizes]

    return labels, sizes


def _percentages(values: Sequence[float]) -> List[float]:
    total = float(sum(values))
    if total == 0.0:
        return [0.0 for _ in values]
    return [v * 100.0 / total for v in values]


def plot_pie(labels: Sequence[str], raw_percents: Sequence[float], colors: dict, labeled: bool, use_log: bool) -> plt.Figure:
    """Create a pie chart figure.

    Parameters
    ----------
    labels : Sequence[str]
        Segment labels in plotting order.
    raw_percents : Sequence[float]
        Segment percentages (0-100), untransformed.
    colors : dict
        Mapping with keys: non_coding, inaccessible_coding, accessible_coding
    labeled : bool
        If True, annotate wedges with raw percentages; else show no wedge text.
    use_log : bool
        If True, wedges are sized by log10(raw_percent + 1).

    Returns
    -------
    matplotlib.figure.Figure
        The created figure.
    """
    fig, ax = plt.subplots(figsize=(5, 5), dpi=150)

    color_list = [
        colors["non_coding_inaccessible"],
        colors["non_coding_accessible"],
        colors["coding_inaccessible"],
        colors["coding_accessible"],
    ]

    # Sizes for plotting (optionally log-transformed)
    sizes_for_plot = [math.log10(p + 1.0) for p in raw_percents] if use_log else list(raw_percents)

    autopct = "%1.1f%%" if labeled else None
    pie_result = ax.pie(
        sizes_for_plot,
        labels=None,  # use legend for labels
        colors=color_list,
        startangle=90,
        counterclock=False,
        autopct=autopct,
        pctdistance=0.7,
    )
    wedges = pie_result[0]
    # If labeled, override the autotexts to show RAW percents, not plot-scaled
    if labeled and len(pie_result) >= 3:
        autotexts = pie_result[2]
        for i, at in enumerate(autotexts):
            if i < len(raw_percents):
                at.set_text(f"{raw_percents[i]:.1f}%")

    # Legend uses raw percentages, not plot-scaled values
    legend_labels = [f"{lab} ({p:.1f}%)" for lab, p in zip(labels, raw_percents)]
    ax.legend(wedges, legend_labels, loc="center left", bbox_to_anchor=(1, 0.5))

    ax.set_aspect("equal")
    if use_log:
        fig.text(
            0.5,
            0.02,
            "Note: slice sizes use log10 transform; percentages reflect raw values",
            ha="center",
            va="center",
            fontsize=8,
        )
    plt.tight_layout()
    return fig


def save_figure(fig: plt.Figure, base_path: Path) -> None:
    """Save figure as PNG and SVG with a shared base path (without suffix).

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.
    base_path : pathlib.Path
        Path without suffix; .png and .svg will be appended.
    """
    png_path = base_path.with_suffix(".png")
    svg_path = base_path.with_suffix(".svg")
    fig.savefig(png_path, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")


def ensure_output_dir(path: str | Path) -> Path:
    """Ensure output directory exists and return it as a Path."""
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Plot genome accessibility pie charts from a summary CSV. "
            "Wedge sizes optionally use log10 transform while labels show raw percentages."
        )
    )
    parser.add_argument("input_csv", help="Path to input CSV file")
    parser.add_argument(
        "output_dir",
        help="Base output directory; figures will be written under <output_dir>/figures",
    )
    parser.add_argument(
        "--rows",
        dest="rows",
        default=None,
        help="Comma-separated list of Display Name rows to include; default uses all rows",
    )
    parser.add_argument(
        "--no-log",
        dest="no_log",
        action="store_true",
        help="Disable log10 sizing of wedges (use raw percentages for sizes)",
    )
    return parser.parse_args()


def main() -> None:
    """Entry point: load data, compute, and write figures for selected rows."""
    args = parse_args()

    input_csv = args.input_csv
    base_output_dir = ensure_output_dir(args.output_dir)
    output_dir = ensure_output_dir(base_output_dir / "figures")

    # Determine selection and log behavior
    if args.rows is None or args.rows.strip() == "":
        display_selection = None
    else:
        display_selection = [r.strip() for r in args.rows.split(",") if r.strip()]
    use_log = not args.no_log

    df = read_summary_csv(input_csv)

    rows = dataframe_to_rows(df, display_selection)

    for row in rows:
        labels, raw_percents = compute_segments(row)

        # Labeled pie
        fig_labeled = plot_pie(labels, raw_percents, COLORS, labeled=True, use_log=use_log)
        base_labeled = output_dir / f"{sanitize_display_name(row.display_name)}_accessibility_pie"
        save_figure(fig_labeled, base_labeled)
        plt.close(fig_labeled)

        # Clear pie (legend only)
        fig_clear = plot_pie(labels, raw_percents, COLORS, labeled=False, use_log=use_log)
        base_clear = output_dir / f"{sanitize_display_name(row.display_name)}_accessibility_pie_clear"
        save_figure(fig_clear, base_clear)
        plt.close(fig_clear)


if __name__ == "__main__":
    main()


