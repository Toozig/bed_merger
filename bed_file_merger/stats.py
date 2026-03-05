from __future__ import annotations

"""Compute BED segment statistics using pandas/numpy."""

from typing import Dict

import numpy as np
import pandas as pd


def compute_basic_stats(bed_df: pd.DataFrame) -> Dict[str, float]:
    if bed_df.empty:
        return {
            "n_peaks": 0,
            "total_bp": 0,
            "min_length": 0,
            "max_length": 0,
            "mean_length": 0.0,
            "median_length": 0.0,
        }

    lengths = (bed_df["end"] - bed_df["start"]).to_numpy()
    return {
        "n_peaks": int(len(lengths)),
        "total_bp": int(np.sum(lengths)),
        "min_length": int(np.min(lengths)),
        "max_length": int(np.max(lengths)),
        "mean_length": float(np.mean(lengths)),
        "median_length": float(np.median(lengths)),
    }


