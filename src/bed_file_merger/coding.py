from __future__ import annotations

"""Coding-region extraction and overlap, using bedtools."""

from pathlib import Path
from typing import Tuple

import os
import subprocess


def check_bedtools_available() -> bool:
    return os.system("command -v bedtools >/dev/null 2>&1") == 0


def extract_coding_regions_from_gtf(gtf_file: Path, output_bed: Path) -> Tuple[int, int]:
    """Extract CDS from GTF to BED (merged), returns (regions, total_bp).

    Uses gtfparse if available to produce a temp BED, then merges with bedtools.
    Otherwise falls back to awk + bedtools merge.
    """
    if not check_bedtools_available():
        raise RuntimeError("bedtools is required for coding-region extraction.")

    try:
        import gtfparse  # type: ignore
        use_basic = False
    except Exception:
        use_basic = True

    output_bed = Path(output_bed)
    output_bed.parent.mkdir(parents=True, exist_ok=True)

    if not use_basic:
        import pandas as pd
        try:
            import polars as pl  # type: ignore
        except Exception:
            pl = None  # type: ignore

        df = gtfparse.read_gtf(str(gtf_file))

        # Handle Polars or pandas
        if pl is not None and hasattr(df, "to_pandas") and str(type(df)).startswith("<class 'polars."):
            cds_pl = df.filter(pl.col("feature") == "CDS")
            if cds_pl.height == 0:
                output_bed.write_text("")
            else:
                bed_pl = cds_pl.select([
                    pl.col("seqname").alias("chr"),
                    (pl.col("start") - 1).alias("start"),
                    pl.col("end").alias("end"),
                ]).sort(["chr", "start"])  # type: ignore
                bed_df = bed_pl.to_pandas()
                tmp_path = str(output_bed) + ".tmp"
                bed_df.to_csv(tmp_path, sep="\t", header=False, index=False)
                with open(output_bed, "w") as out_handle:
                    subprocess.run(["bedtools", "merge", "-i", tmp_path], check=True, stdout=out_handle)
                os.unlink(tmp_path)
        else:
            # Assume pandas
            cds_df = df[df["feature"] == "CDS"].copy()
            if cds_df.empty:
                output_bed.write_text("")
            else:
                bed_df = pd.DataFrame({
                    "chr": cds_df["seqname"],
                    "start": cds_df["start"] - 1,
                    "end": cds_df["end"],
                }).sort_values(["chr", "start"])  # type: ignore
                tmp_path = str(output_bed) + ".tmp"
                bed_df.to_csv(tmp_path, sep="\t", header=False, index=False)
                with open(output_bed, "w") as out_handle:
                    subprocess.run(["bedtools", "merge", "-i", tmp_path], check=True, stdout=out_handle)
                os.unlink(tmp_path)
    else:
        cmd = (
            "awk -F'\\t' '$3 == \"CDS\" {print $1 \"\\t\" ($4-1) \"\\t\" $5}' "
            f"{gtf_file} | sort -k1,1 -k2,2n | bedtools merge -i - > {output_bed}"
        )
        subprocess.run(["bash", "-c", cmd], check=True)

    # stats
    wc = subprocess.run(["wc", "-l", str(output_bed)], capture_output=True, text=True, check=True)
    regions = int(wc.stdout.split()[0])
    awk = subprocess.run(["awk", "-F\t", "{sum += $3 - $2} END {print sum}", str(output_bed)],
                         capture_output=True, text=True, check=True)
    total_bp = int(awk.stdout.strip() or "0")
    return regions, total_bp


def count_intersections_and_bp(a_bed_path: Path, b_bed_path: Path) -> Tuple[int, int]:
    """Return (segments_overlapping, bp_overlapping) using bedtools intersect."""
    if not check_bedtools_available():
        return 0, 0

    # unique segments intersecting
    res = subprocess.run([
        "bedtools", "intersect", "-a", str(a_bed_path), "-b", str(b_bed_path), "-u"
    ], capture_output=True, text=True)
    if res.returncode == 0 and res.stdout.strip():
        segments_overlapping = len([ln for ln in res.stdout.splitlines() if ln.strip()])
    else:
        segments_overlapping = 0

    res_bp = subprocess.run([
        "bedtools", "intersect", "-a", str(a_bed_path), "-b", str(b_bed_path)
    ], capture_output=True, text=True)
    if res_bp.returncode == 0 and res_bp.stdout.strip():
        bp = 0
        for ln in res_bp.stdout.splitlines():
            parts = ln.split("\t")
            if len(parts) >= 3:
                bp += int(parts[2]) - int(parts[1])
        overlap_bp = bp
    else:
        overlap_bp = 0

    return segments_overlapping, overlap_bp


