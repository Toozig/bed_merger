from __future__ import annotations

"""I/O helpers for reading BED files into DataFrames."""

from pathlib import Path
from typing import Tuple, Optional, List

import pandas as pd


def load_bed_to_dataframe(bed_path: Path, extra_column_names: Optional[List[str]] = None) -> pd.DataFrame:
    """Load a BED file into a DataFrame with columns chr, start, end, plus extras if present.

    Assumes tab-delimited; tolerates spaces. Filters malformed rows and ensures integer coords.
    """
    with open(bed_path, "r", newline="") as handle:
        lines = [line.strip() for line in handle if line.strip() and not line.startswith("#")]

    rows = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 3:
            parts = line.split()
        if len(parts) >= 3:
            rows.append(parts)

    if not rows:
        return pd.DataFrame(columns=["chr", "start", "end"])  # empty

    df = pd.DataFrame(rows)
    base_cols = ["chr", "start", "end"]
    if df.shape[1] > 3:
        expected_extras = df.shape[1] - 3
        if extra_column_names is not None:
            if len(extra_column_names) < expected_extras:
                # not enough provided, fill remaining with default names
                provided = list(extra_column_names)
                missing = expected_extras - len(provided)
                extras = provided + [f"col_{i}" for i in range(4 + len(provided), 4 + len(provided) + missing)]
            elif len(extra_column_names) > expected_extras:
                raise ValueError(
                    f"Too many column names provided: got {len(extra_column_names)} for file {bed_path.name}, "
                    f"but file has only {expected_extras} extra columns"
                )
            else:
                extras = list(extra_column_names)
        else:
            extras = [f"col_{i}" for i in range(4, df.shape[1] + 1)]
        columns = base_cols + extras
    else:
        columns = base_cols[: df.shape[1]]
    df.columns = columns

    if "start" in df.columns and "end" in df.columns:
        df["start"] = pd.to_numeric(df["start"], errors="coerce")
        df["end"] = pd.to_numeric(df["end"], errors="coerce")
        df = df.dropna(subset=["start", "end"]).copy()
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)
        df = df[df["end"] > df["start"]]

    return df.reset_index(drop=True)


