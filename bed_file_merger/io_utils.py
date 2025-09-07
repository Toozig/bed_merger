from __future__ import annotations

"""I/O helpers for reading BED files into DataFrames."""

from pathlib import Path
from typing import Tuple, Optional, List

import pandas as pd


def read_bed_file(bed_path: Path) -> pd.DataFrame:
    """Read a BED file using pandas with tab delimiter, ignoring comments and blanks.

    Raises ValueError if the file contains no data rows.
    """
    df = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        comment="#",
        skip_blank_lines=True,
        dtype=str,
        engine="python",
    )
    # Drop fully empty rows (in case of stray separators)
    df = df.dropna(how="all")
    if df.shape[0] == 0:
        raise ValueError(f"BED file '{bed_path}' contains no data rows after filtering comments/blanks.")
    return df


def normalize_bed_columns(df: pd.DataFrame, extra_column_names: Optional[List[str]] = None) -> pd.DataFrame:
    """Normalize a raw BED DataFrame to have chr,start,end plus optional extras.

    - Requires at least 3 columns; otherwise raises ValueError.
    - Renames first three columns to chr,start,end.
    - Validates and assigns extra column names.
    - Converts start/end to integers and filters invalid rows (end > start).
    """
    if df.shape[1] < 3:
        raise ValueError("BED file has fewer than 3 columns (chr, start, end) after parsing.")

    base_cols = ["chr", "start", "end"]
    n_extra = max(df.shape[1] - 3, 0)
    if n_extra > 0:
        if extra_column_names is not None:
            if len(extra_column_names) > n_extra:
                raise ValueError(
                    f"Too many extra column names provided: {len(extra_column_names)} but file has only {n_extra} extra columns"
                )
            provided = list(extra_column_names)
            missing = n_extra - len(provided)
            extras = provided + [f"col_{i}" for i in range(4 + len(provided), 4 + len(provided) + missing)]
        else:
            extras = [f"col_{i}" for i in range(4, 4 + n_extra)]
        columns = base_cols + extras
    else:
        columns = base_cols[: df.shape[1]]

    # Assign column names
    df = df.copy()
    df.columns = columns

    # Coerce start/end to numeric and clean
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"]).copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df = df[df["end"] > df["start"]]
    return df.reset_index(drop=True)


def load_bed_to_dataframe(bed_path: Path, extra_column_names: Optional[List[str]] = None) -> pd.DataFrame:
    """Public API: read a BED and normalize its columns."""
    raw = read_bed_file(bed_path)
    return normalize_bed_columns(raw, extra_column_names=extra_column_names)


