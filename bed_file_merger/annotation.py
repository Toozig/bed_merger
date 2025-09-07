from __future__ import annotations

"""Peak annotation utilities built on top of bedtools.

Creates a BED with genome annotations from a GTF and annotates peaks by
intersecting with the annotation BED using a priority order.
"""

from pathlib import Path
from typing import List, Dict

import os
import subprocess
import pandas as pd


def check_bedtools_available() -> bool:
    return os.system("command -v bedtools >/dev/null 2>&1") == 0


def build_annotation_bed_from_gtf(gtf_file: Path, out_bed: Path) -> Path:
    """Build a BED-like file: chr, start(0-based), end, feature from a GTF.

    Uses gtfparse when available; falls back to awk. The BED is sorted.
    """
    out_bed = Path(out_bed)
    out_bed.parent.mkdir(parents=True, exist_ok=True)

    try:
        import gtfparse  # type: ignore
        df = gtfparse.read_gtf(str(gtf_file))
        if df.shape[0] == 0:
            out_bed.write_text("")
            return out_bed
        # Produce required columns
        ann_df = pd.DataFrame({
            "chr": df["seqname"],
            "start": df["start"] - 1,
            "end": df["end"],
            "feature": df["feature"],
        })
        ann_df = ann_df.sort_values(["chr", "start"])  # type: ignore
        ann_df.to_csv(out_bed, sep="\t", header=False, index=False)
        return out_bed
    except Exception:
        # Fallback awk: chr, start-1, end, feature
        cmd = (
            f"awk -F'\\t' '{{print $1 \"\\t\" ($4-1) \"\\t\" $5 \"\\t\" $3}}' {gtf_file} | "
            f"sort -k1,1 -k2,2n > {out_bed}"
        )
        subprocess.run(["bash", "-c", cmd], check=True)
        return out_bed


def annotate_peaks_with_gtf(
    peaks_df: pd.DataFrame,
    annotation_bed: Path,
    priority: List[str],
    tmp_dir: Path,
) -> List[str]:
    """Annotate each peak by intersecting with the annotation BED and picking the
    highest priority feature. Returns a list of annotations ('' if none).
    """
    if peaks_df.empty:
        return [""] * 0

    if not check_bedtools_available():
        return [""] * len(peaks_df)

    # Create a temporary bed with row ids to map back
    tmp_dir = Path(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    tmp_bed = tmp_dir / "_peaks_for_annot.tmp.bed"
    peaks_with_id = peaks_df[["chr", "start", "end"]].copy()
    peaks_with_id["row_id"] = range(len(peaks_with_id))
    peaks_with_id.to_csv(tmp_bed, sep="\t", header=False, index=False)

    # intersect -wa -wb to get both A (peaks) and B (annotation)
    cmd = [
        "bedtools", "intersect",
        "-a", str(tmp_bed),
        "-b", str(annotation_bed),
        "-wa", "-wb",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    annotations: Dict[int, List[str]] = {}
    if result.returncode == 0 and result.stdout.strip():
        for line in result.stdout.splitlines():
            # A: chr, start, end, row_id  B: chr, start, end, feature
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            try:
                row_id = int(parts[3])
            except ValueError:
                continue
            feature = parts[7]
            if row_id not in annotations:
                annotations[row_id] = []
            annotations[row_id].append(feature)

    try:
        tmp_bed.unlink(missing_ok=True)
    except Exception:
        pass

    # Priority mapping (lower index => higher priority)
    priority_index = {name: idx for idx, name in enumerate(priority)}

    picked: List[str] = []
    for idx in range(len(peaks_df)):
        feats = annotations.get(idx, [])
        if not feats:
            picked.append("")
            continue
        # keep only features that are in the priority list
        filtered = [f for f in feats if f in priority_index]
        if not filtered:
            # only non-priority annotations (e.g., transcript) -> treat as no annotation
            picked.append("")
            continue
        # choose min index among filtered
        best = min(filtered, key=lambda f: priority_index[f])
        picked.append(best)

    return picked


def derive_coding_bed_from_annotation(annotation_bed: Path, output_bed: Path) -> Path:
    """Create a coding BED (merged CDS) from an annotation BED built from the GTF.

    This avoids re-reading the GTF by filtering the annotation BED for CDS and
    merging intervals with bedtools.
    """
    output_bed = Path(output_bed)
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    # Filter for CDS (4th column is feature) and merge
    cmd = (
        f"awk -F'\t' '$4 == \"CDS\" {{print $1 \"\t\" $2 \"\t\" $3}}' {annotation_bed} | "
        f"sort -k1,1 -k2,2n | bedtools merge -i - > {output_bed}"
    )
    subprocess.run(["bash", "-c", cmd], check=True)
    return output_bed


