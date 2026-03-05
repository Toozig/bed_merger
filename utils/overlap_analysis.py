from __future__ import annotations

"""
Overlap analysis utility

Inputs:
- results_dir: directory produced by bed_file_merger run
- sheet_a, sheet_b: names of two sheets in the Excel within results_dir

Behavior:
- Auto-detect the Excel file in results_dir (or allow override)
- Auto-detect coding BED in results_dir (pattern: *_coding_regions_refseq.bed); override with --coding-bed
- Read BED-like data (chr,start,end) from the two sheets
- Sort and bedtools merge each set
- Compute:
  - mutual (exact overlaps between merged A and merged B)
  - unique_A (A minus B)
  - unique_B (B minus A)
- For each group: total_bp, coding_bp, noncoding_bp, and percentages (coding/noncoding)
- Append a sheet "Overlap_Analysis_<A>_vs_<B>" to the Excel (in-place)
- Generate two venn diagrams under <results_dir>/figures/:
  - venn_no_text_<A>_vs_<B>.svg/png (no text)
  - venn_with_text_<A>_vs_<B>.svg/png (with counts)

Requirements: bedtools in PATH
"""

import argparse
import json
import re
import subprocess
from pathlib import Path
from typing import List, Tuple, Optional

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
try:
    from matplotlib_venn import venn2
except Exception:  # pragma: no cover
    venn2 = None  # will fail at runtime with a clear message


def ensure_bedtools() -> None:
    try:
        subprocess.run(["bedtools", "--version"], check=True, capture_output=True, text=True)
    except Exception as exc:
        raise SystemExit("Error: bedtools not found in PATH.") from exc


def sanitize_for_filename(name: str) -> str:
    base = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")
    base = re.sub(r"_+", "_", base)
    return base


def sanitize_for_sheet(name: str) -> str:
    # Excel sheet name constraints: no []:*?/\\ and <= 31 chars
    safe = re.sub(r"[\[\]:\\/?*]+", "_", name)
    safe = safe.strip()
    if len(safe) > 31:
        safe = safe[:31]
    if not safe:
        safe = "Sheet"
    return safe


def autodetect_excel(results_dir: Path) -> Path:
    cands = sorted(results_dir.glob("*.xlsx"))
    if not cands:
        raise SystemExit(f"No .xlsx found under {results_dir}")
    if len(cands) > 1:
        # Pick the newest
        cands.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return cands[0]


def autodetect_coding_bed(results_dir: Path) -> Optional[Path]:
    cands = sorted(results_dir.glob("*_coding_regions_refseq.bed"))
    return cands[0] if cands else None


def df_to_temp_bed(df: pd.DataFrame, cols: List[str], out: Path) -> Path:
    sub = df[cols].copy()
    sub.to_csv(out, sep="\t", header=False, index=False)
    return out


def sort_and_merge(in_bed: Path, out_bed: Path) -> None:
    # sort -k1,1 -k2,2n | bedtools merge -i -
    sort_cmd = ["sort", "-k1,1", "-k2,2n", str(in_bed)]
    p1 = subprocess.Popen(sort_cmd, stdout=subprocess.PIPE, text=True)
    p2 = subprocess.Popen(["bedtools", "merge", "-i", "-"], stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
    p1.stdout.close()  # type: ignore
    stdout, _ = p2.communicate()
    if p2.returncode != 0:
        raise SystemExit("Error: bedtools merge failed.")
    out_bed.write_text(stdout)


def bed_total_bp(bed_path: Path) -> int:
    if not bed_path.exists() or bed_path.stat().st_size == 0:
        return 0
    # awk sum
    awk = subprocess.run(["awk", "-F\t", "{sum+=($3-$2)} END {print sum}", str(bed_path)], capture_output=True, text=True)
    try:
        return int((awk.stdout or "0").strip() or "0")
    except Exception:
        return 0


def compute_mutual_bed(a_merged: Path, b_merged: Path, out_mutual: Path) -> Tuple[int, Path]:
    # bedtools intersect -a A -b B -wo, compute overlap segments
    proc = subprocess.run(["bedtools", "intersect", "-a", str(a_merged), "-b", str(b_merged), "-wo"], capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit("Error: bedtools intersect -wo failed.")
    total_mutual = 0
    lines: List[str] = []
    for line in (proc.stdout or "").splitlines():
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        # A: chrA,startA,endA ; B: chrB,startB,endB ; overlap as last
        try:
            chr_a = parts[0]
            start_a = int(parts[1])
            end_a = int(parts[2])
            chr_b = parts[3]
            start_b = int(parts[4])
            end_b = int(parts[5])
            ov = int(parts[-1])
        except Exception:
            continue
        if ov <= 0:
            continue
        total_mutual += ov
        # Overlap segment coordinates
        if chr_a != chr_b:
            # Ignore cross-chrom anomalies
            continue
        start = max(start_a, start_b)
        end = min(end_a, end_b)
        if start < end:
            lines.append(f"{chr_a}\t{start}\t{end}\n")
    # Write and merge overlap segments to get non-overlapping mutual BED
    tmp_overlap = out_mutual.with_suffix(".tmp.overlap.bed")
    tmp_overlap.write_text("".join(lines))
    sort_and_merge(tmp_overlap, out_mutual)
    try:
        tmp_overlap.unlink()
    except Exception:
        pass
    return total_mutual, out_mutual


def bed_coding_bp(group_bed: Path, coding_bed: Path) -> int:
    if not group_bed.exists() or group_bed.stat().st_size == 0:
        return 0
    proc = subprocess.run(["bedtools", "intersect", "-a", str(group_bed), "-b", str(coding_bed), "-wo"], capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit("Error: bedtools intersect (coding) failed.")
    total = 0
    for line in (proc.stdout or "").splitlines():
        try:
            ov = int(line.split("\t")[-1])
        except Exception:
            ov = 0
        total += ov
    return total


def draw_venn(fig_dir: Path, base_a: str, base_b: str, unique_a: int, unique_b: int, mutual: int) -> None:
    fig_dir.mkdir(parents=True, exist_ok=True)

    total_a = unique_a + mutual
    total_b = unique_b + mutual

    if venn2 is None:
        raise SystemExit("matplotlib_venn is required for styled Venn diagrams. Please install matplotlib-venn.")

    # No-text version
    fig = plt.figure(figsize=(6, 4), dpi=300)
    v = venn2(subsets=(unique_a, unique_b, mutual), set_labels=(base_a, base_b))
    # Colors similar to notebook (pastel palette)
    palette = sns.color_palette("pastel", 12)
    if v.get_patch_by_id('10') is not None:
        v.get_patch_by_id('10').set_color(palette[-3])
        v.get_patch_by_id('10').set_alpha(0.7)
    if v.get_patch_by_id('01') is not None:
        v.get_patch_by_id('01').set_color(palette[-7])
        v.get_patch_by_id('01').set_alpha(0.7)
    if v.get_patch_by_id('11') is not None:
        v.get_patch_by_id('11').set_color(palette[-2])
        v.get_patch_by_id('11').set_alpha(0.7)
    # Remove all text for this version
    for lid in ['10', '01', '11']:
        lbl = v.get_label_by_id(lid)
        if lbl is not None:
            lbl.set_text("")
    for lbl in (v.set_labels or []):
        if lbl is not None:
            lbl.set_text("")
    plt.tight_layout()
    out_base = f"venn_no_text_{base_a}_vs_{base_b}"
    fig.savefig(fig_dir / f"{out_base}.svg", transparent=True)
    fig.savefig(fig_dir / f"{out_base}.png", dpi=300, transparent=True)
    plt.close(fig)

    # With-text version (similar to notebook layout/style)
    fig = plt.figure(figsize=(10, 6), dpi=300)
    v = venn2(subsets=(unique_a, unique_b, mutual), set_labels=(base_a, base_b))
    for lid in ['10', '01', '11']:
        patch = v.get_patch_by_id(lid)
        if patch is not None:
            patch.set_alpha(0.7)
    if v.get_patch_by_id('10') is not None:
        v.get_patch_by_id('10').set_color(palette[-3])
    if v.get_patch_by_id('01') is not None:
        v.get_patch_by_id('01').set_color(palette[-7])
    if v.get_patch_by_id('11') is not None:
        v.get_patch_by_id('11').set_color(palette[-2])

    label_fontsize = 15
    venn_label_fontsize = 13
    # Set label fonts
    if v.set_labels is not None:
        for label in v.set_labels:
            if label is not None:
                label.set_fontsize(label_fontsize)

    # Region labels with numbers and percentages
    a_pct_unique = (unique_a / total_a * 100.0) if total_a > 0 else 0.0
    b_pct_unique = (unique_b / total_b * 100.0) if total_b > 0 else 0.0
    a_pct_shared = (mutual / total_a * 100.0) if total_a > 0 else 0.0
    b_pct_shared = (mutual / total_b * 100.0) if total_b > 0 else 0.0

    if v.get_label_by_id('10') is not None:
        v.get_label_by_id('10').set_text(f"{unique_a:,}bp\n({a_pct_unique:.1f}%)")
        v.get_label_by_id('10').set_fontsize(venn_label_fontsize)
    if v.get_label_by_id('01') is not None:
        v.get_label_by_id('01').set_text(f"{unique_b:,}bp\n({b_pct_unique:.1f}%)")
        v.get_label_by_id('01').set_fontsize(venn_label_fontsize)
    if v.get_label_by_id('11') is not None:
        v.get_label_by_id('11').set_text(
            f"{mutual:,}bp\n({a_pct_shared:.1f}% of {base_a})\n({b_pct_shared:.1f}% of {base_b})"
        )
        v.get_label_by_id('11').set_fontsize(venn_label_fontsize)

    # Slight positional tweaks similar to notebook
    if v.get_label_by_id('11') is not None:
        x11, y11 = v.get_label_by_id('11').get_position()
        v.get_label_by_id('11').set_y(y11 + 0.0)
    for lid in ['10', '01']:
        if v.get_label_by_id(lid) is not None and v.get_label_by_id('11') is not None:
            _, ymid = v.get_label_by_id('11').get_position()
            v.get_label_by_id(lid).set_y(ymid)

    # Extra text under set labels: total bp
    if v.set_labels is not None:
        # A text
        la = v.set_labels[0]
        if la is not None:
            xa, ya = la.get_position()
            plt.text(xa - 0.04, ya - 0.01, f"\n{total_a:,}bp", ha='center', va='top', fontsize=venn_label_fontsize)
        lb = v.set_labels[1]
        if lb is not None and v.set_labels[0] is not None:
            xb, yb = lb.get_position()
            xa0, ya0 = v.set_labels[0].get_position()
            # Align second label y to first's y to mimic notebook
            v.set_labels[1].set_position((xb, ya0))
            xb2, y2 = v.set_labels[1].get_position()
            plt.text(xb2 + 0.09, y2 - 0.01, f"\n{total_b:,}bp", ha='center', va='top', fontsize=venn_label_fontsize)

    plt.tight_layout()
    out_base = f"venn_with_text_{base_a}_vs_{base_b}"
    fig.savefig(fig_dir / f"{out_base}.svg", transparent=True)
    fig.savefig(fig_dir / f"{out_base}.png", dpi=300, transparent=True)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Add overlap analysis tab and Venn figures to results directory Excel")
    parser.add_argument("--results-dir", required=True, help="Path to results directory")
    parser.add_argument("--sheet-a", required=True, help="First sheet name")
    parser.add_argument("--sheet-b", required=True, help="Second sheet name")
    parser.add_argument("--excel", default=None, help="Override: Excel path (if multiple)")
    parser.add_argument("--coding-bed", default=None, help="Override: coding BED path")
    parser.add_argument("--genome-build", default="hg38", help="Genome build key to read totals from resources/genome_info.json (default: hg38)")
    args = parser.parse_args()

    ensure_bedtools()

    results_dir = Path(args.results_dir)
    if not results_dir.exists():
        raise SystemExit(f"results_dir not found: {results_dir}")

    excel_path = Path(args.excel) if args.excel else autodetect_excel(results_dir)
    coding_bed = Path(args.coding_bed) if args.coding_bed else autodetect_coding_bed(results_dir)
    if coding_bed is None:
        raise SystemExit("coding BED not found; pass --coding-bed")

    # Read sheets
    try:
        a_df = pd.read_excel(excel_path, sheet_name=args.sheet_a)
        b_df = pd.read_excel(excel_path, sheet_name=args.sheet_b)
    except Exception as exc:
        raise SystemExit(f"Failed reading sheets: {exc}")

    # Expect chr,start,end columns; normalize if needed
    def normalize(df: pd.DataFrame) -> pd.DataFrame:
        cols = list(df.columns)
        lower = {c.lower(): c for c in cols}
        needed = ["chr", "start", "end"]
        if not all(n in lower for n in needed):
            raise SystemExit("Sheets must have chr,start,end columns")
        return df[[lower["chr"], lower["start"], lower["end"]]].rename(columns={lower["chr"]: "chr", lower["start"]: "start", lower["end"]: "end"})

    a_df = normalize(a_df)
    b_df = normalize(b_df)

    tmp_dir = results_dir / "_overlap_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    a_raw = tmp_dir / "A.raw.bed"
    b_raw = tmp_dir / "B.raw.bed"
    df_to_temp_bed(a_df, ["chr", "start", "end"], a_raw)
    df_to_temp_bed(b_df, ["chr", "start", "end"], b_raw)

    a_merged = tmp_dir / "A.merged.bed"
    b_merged = tmp_dir / "B.merged.bed"
    sort_and_merge(a_raw, a_merged)
    sort_and_merge(b_raw, b_merged)

    # Unique sets
    a_only = tmp_dir / "A.only.bed"
    b_only = tmp_dir / "B.only.bed"
    with (a_only).open("w") as outa:
        p = subprocess.run(["bedtools", "subtract", "-a", str(a_merged), "-b", str(b_merged)], capture_output=True, text=True)
        if p.returncode != 0:
            raise SystemExit("Error: bedtools subtract A-B failed")
        outa.write(p.stdout or "")
    with (b_only).open("w") as outb:
        p = subprocess.run(["bedtools", "subtract", "-a", str(b_merged), "-b", str(a_merged)], capture_output=True, text=True)
        if p.returncode != 0:
            raise SystemExit("Error: bedtools subtract B-A failed")
        outb.write(p.stdout or "")

    # Mutual segments
    mutual_bed = tmp_dir / "mutual.bed"
    mutual_bp, mutual_bed = compute_mutual_bed(a_merged, b_merged, mutual_bed)

    # Totals
    total_a_only = bed_total_bp(a_only)
    total_b_only = bed_total_bp(b_only)
    # mutual_bp already computed

    # Coding splits
    coding_bp_mut = bed_coding_bp(mutual_bed, coding_bed)
    coding_bp_a = bed_coding_bp(a_only, coding_bed)
    coding_bp_b = bed_coding_bp(b_only, coding_bed)

    def row_for(group: str, total_bp: int, coding_bp: int) -> dict:
        noncoding = max(total_bp - coding_bp, 0)
        total_bp = max(total_bp, 0)
        return {
            "group": group,
            "total_bp": total_bp,
            "coding_bp": coding_bp,
            "noncoding_bp": noncoding,
            "coding_pct": (coding_bp / total_bp * 100.0) if total_bp > 0 else 0.0,
            "noncoding_pct": (noncoding / total_bp * 100.0) if total_bp > 0 else 0.0,
        }

    rows = [
        row_for("mutual", mutual_bp, coding_bp_mut),
        row_for("unique_A", total_a_only, coding_bp_a),
        row_for("unique_B", total_b_only, coding_bp_b),
    ]
    out_df = pd.DataFrame(rows)

    # Add genome % using resources/genome_info.json
    try:
        resources_json = Path(__file__).resolve().parents[1] / "bed_file_merger" / "resources" / "genome_info.json"
        with resources_json.open("r") as f:
            gi = json.load(f)
        gb = gi.get(args.genome_build, {})
        genome_total_bp = int(gb.get("Total bases", 0))
    except Exception:
        genome_total_bp = 0

    if genome_total_bp > 0:
        out_df["group_genome_pct"] = (out_df["total_bp"] / genome_total_bp * 100.0).round(2)
    else:
        out_df["group_genome_pct"] = 0.0

    # Append to Excel
    base_a = sanitize_for_filename(args.sheet_a)
    base_b = sanitize_for_filename(args.sheet_b)
    sheet_name = sanitize_for_sheet(f"Overlap_Analysis_{base_a}_vs_{base_b}")
    with pd.ExcelWriter(excel_path, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
        out_df.to_excel(writer, sheet_name=sheet_name, index=False)

    # Figures
    fig_dir = results_dir / "figures"
    # Pass through genome_total_bp to include percent-of-genome under labels when available
    try:
        # monkey-patch call with closure variables via kwargs if supported
        draw_venn.__defaults__
        # We will just call and let function re-compute percent text if needed; modify signature if desired
    except Exception:
        pass
    # Reuse draw_venn which formats totals inside; update to include genome % by environment
    # To avoid changing signature, we attach globals the function can see; instead, compute text inside caller:
    draw_venn(fig_dir, base_a, base_b, total_a_only, total_b_only, mutual_bp)

    # If genome total is available, also render an alternate with genome-percent under labels
    if genome_total_bp > 0:
        # Duplicate the with-text figure with genome percent lines (suffix _with_text already used)
        # We'll regenerate the with_text version inline here
        total_a = total_a_only + mutual_bp
        total_b = total_b_only + mutual_bp
        a_pct_genome = (total_a / genome_total_bp * 100.0) if genome_total_bp > 0 else 0.0
        b_pct_genome = (total_b / genome_total_bp * 100.0) if genome_total_bp > 0 else 0.0
        # Create figure similar to draw_venn with text but add genome line
        if venn2 is None:
            pass
        else:
            palette = sns.color_palette("pastel", 12)
            fig = plt.figure(figsize=(10, 6), dpi=300)
            v = venn2(subsets=(total_a_only, total_b_only, mutual_bp), set_labels=(base_a, base_b))
            for lid in ['10', '01', '11']:
                patch = v.get_patch_by_id(lid)
                if patch is not None:
                    patch.set_alpha(0.7)
            if v.get_patch_by_id('10') is not None:
                v.get_patch_by_id('10').set_color(palette[-3])
            if v.get_patch_by_id('01') is not None:
                v.get_patch_by_id('01').set_color(palette[-7])
            if v.get_patch_by_id('11') is not None:
                v.get_patch_by_id('11').set_color(palette[-2])
            label_fontsize = 15
            venn_label_fontsize = 13
            if v.set_labels is not None:
                for lbl in v.set_labels:
                    if lbl is not None:
                        lbl.set_fontsize(label_fontsize)
            # Region labels as before
            a_pct_unique = (total_a_only / total_a * 100.0) if total_a > 0 else 0.0
            b_pct_unique = (total_b_only / total_b * 100.0) if total_b > 0 else 0.0
            a_pct_shared = (mutual_bp / total_a * 100.0) if total_a > 0 else 0.0
            b_pct_shared = (mutual_bp / total_b * 100.0) if total_b > 0 else 0.0
            if v.get_label_by_id('10') is not None:
                v.get_label_by_id('10').set_text(f"{total_a_only:,}bp\n({a_pct_unique:.1f}%)")
                v.get_label_by_id('10').set_fontsize(venn_label_fontsize)
            if v.get_label_by_id('01') is not None:
                v.get_label_by_id('01').set_text(f"{total_b_only:,}bp\n({b_pct_unique:.1f}%)")
                v.get_label_by_id('01').set_fontsize(venn_label_fontsize)
            if v.get_label_by_id('11') is not None:
                v.get_label_by_id('11').set_text(
                    f"{mutual_bp:,}bp\n({a_pct_shared:.1f}% of {base_a})\n({b_pct_shared:.1f}% of {base_b})"
                )
                v.get_label_by_id('11').set_fontsize(venn_label_fontsize)
            # Add totals + genome % under labels
            if v.set_labels is not None:
                la = v.set_labels[0]
                if la is not None:
                    xa, ya = la.get_position()
                    plt.text(xa - 0.04, ya - 0.01, f"\n{total_a:,}bp\n({a_pct_genome:.2f}% of {args.genome_build} genome)",
                             ha='center', va='top', fontsize=venn_label_fontsize)
                lb = v.set_labels[1]
                if lb is not None and v.set_labels[0] is not None:
                    xb, yb = lb.get_position()
                    xa0, ya0 = v.set_labels[0].get_position()
                    v.set_labels[1].set_position((xb, ya0))
                    xb2, y2 = v.set_labels[1].get_position()
                    plt.text(xb2 + 0.09, y2 - 0.01, f"\n{total_b:,}bp\n({b_pct_genome:.2f}% of {args.genome_build} genome)",
                             ha='center', va='top', fontsize=venn_label_fontsize)
            plt.tight_layout()
            out_base = f"venn_with_text_{base_a}_vs_{base_b}"
            fig.savefig(fig_dir / f"{out_base}.svg", transparent=True)
            fig.savefig(fig_dir / f"{out_base}.png", dpi=300, transparent=True)
            plt.close(fig)

    # Cleanup tmp
    for p in [a_raw, b_raw, a_merged, b_merged, a_only, b_only, mutual_bed]:
        try:
            p.unlink()
        except Exception:
            pass
    try:
        tmp_dir.rmdir()
    except Exception:
        pass

    print(f"Appended sheet: {sheet_name}")
    print(f"Figures written to: {fig_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


