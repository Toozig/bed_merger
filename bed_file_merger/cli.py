from __future__ import annotations

"""CLI for BED file analysis and coding-region operations."""

from pathlib import Path
from typing import List, Optional
import subprocess as sp
import sys
import typer
import pandas as pd
import json
from importlib import resources
from .io_utils import load_bed_to_dataframe, read_bed_file, normalize_bed_columns
from .stats import compute_basic_stats
from .coding import (
    extract_coding_regions_from_gtf,
    count_intersections_and_bp,
    check_bedtools_available,
)
from .report import process_per_file_frames, compute_summary_df, save_excel_report
from .models import AnalysisConfig, BedFileStats, RunConfig
from .annotation import build_annotation_bed_from_gtf, annotate_peaks_with_gtf, derive_coding_bed_from_annotation
import yaml
import shutil
import filecmp


app = typer.Typer(add_completion=False, help="BED analysis CLI")


def _generate_report(
    bed_files: List[Path],
    output_excel: Path,
    genome_build: str,
    refseq_gtf: Optional[Path] = None,
    coding_bed: Optional[Path] = None,
    extra_column_names: Optional[List[Optional[List[str]]]] = None,
    file_names: Optional[List[str]] = None,
    genome_total_bp: Optional[int] = None,
    coding_total_bp: Optional[int] = None,
    tmp_dir: Optional[Path] = None,
    add_id: bool = False,
    id_column_name: str = "id",
    merged_bed_path: Optional[Path] = None,
    merged_sheet_name: str = "merged",
) -> None:
    """Core report generation used by both CLI modes."""
    # Load BEDs
    per_file_dfs: List[pd.DataFrame] = []
    per_file_names: List[str] = []
    extra_names: List[Optional[List[str]]] = [None] * len(bed_files)
    if extra_column_names is not None:
        if len(extra_column_names) != len(bed_files):
            raise typer.Exit(code=1)
        extra_names = extra_column_names

    for idx, p in enumerate(bed_files):
        df = load_bed_to_dataframe(p, extra_column_names=extra_names[idx])
        # Optionally add ID column using provided file_names (skip if already present)
        if add_id and (id_column_name not in df.columns):
            if file_names is None or idx >= len(file_names):
                typer.echo("Error: add_id is enabled but file_names are missing or mismatched.", err=True)
                raise typer.Exit(code=1)
            display_name = file_names[idx]
            df = df.copy()
            df[id_column_name] = [f"{display_name}.{i+1}" for i in range(len(df))]
        per_file_dfs.append(df)
        per_file_names.append(p.name)

    # Optional coding: prepare coding_bed if GTF provided
    coding_bed_path: Optional[Path] = None
    if refseq_gtf is not None:
        if not check_bedtools_available():
            typer.echo("Error: bedtools not found in PATH but required for coding analysis.", err=True)
            raise typer.Exit(code=1)
        # Build one annotation BED once
        ann_bed = Path(output_excel).with_suffix("").parent / "annotation_from_gtf.bed"
        if not ann_bed.exists() or ann_bed.stat().st_size == 0:
            build_annotation_bed_from_gtf(refseq_gtf, ann_bed)
        # Derive coding bed from annotation (avoids re-reading GTF)
        if coding_bed is not None:
            coding_bed_path = coding_bed
            coding_bed_path.parent.mkdir(parents=True, exist_ok=True)
        else:
            coding_bed_path = Path(str(output_excel).replace(".xlsx", "_coding_regions.bed"))
        if not coding_bed_path.exists() or coding_bed_path.stat().st_size == 0:
            derive_coding_bed_from_annotation(ann_bed, coding_bed_path)

    # Helper functions to avoid duplication
    stats_records: List[BedFileStats] = []
    tmp_paths: List[Path] = []

    def _annotate_df(df: pd.DataFrame) -> pd.DataFrame:
        if refseq_gtf is None:
            return df
        ann_bed_local = Path(output_excel).with_suffix("").parent / "annotation_from_gtf.bed"
        priority_local = [
            "3UTR",
            "5UTR",
            "CDS",
        ]
        ann_local = annotate_peaks_with_gtf(
            df,
            ann_bed_local,
            priority_local,
            Path(tmp_dir) if tmp_dir is not None else Path(typer.get_app_dir("bed-file-merger")),
        )
        if len(ann_local) == len(df):
            result = df.copy()
            result["genomic_annotation"] = ann_local
            return result
        return df

    def _compute_stats(df: pd.DataFrame, name_for_tmp: Path, file_name_field: str) -> BedFileStats:
        basic_local = compute_basic_stats(df)
        coding_peaks_local = None
        coding_bp_local = None
        coding_percent_local = None
        if coding_bed_path is not None and not df.empty:
            tmp_base_local = Path(tmp_dir) if tmp_dir is not None else Path(typer.get_app_dir("bed-file-merger"))
            tmp_base_local.mkdir(parents=True, exist_ok=True)
            tmp_local = tmp_base_local / f"tmp_{Path(name_for_tmp).stem}.bed"
            df[["chr", "start", "end"]].to_csv(tmp_local, sep="\t", header=False, index=False)
            tmp_paths.append(tmp_local)
            segs_local, bp_local = count_intersections_and_bp(tmp_local, coding_bed_path)
            coding_peaks_local = segs_local
            coding_bp_local = bp_local
            total_bp_local = basic_local["total_bp"]
            coding_percent_local = (bp_local / total_bp_local * 100.0) if total_bp_local > 0 else 0.0
        return BedFileStats(
            file_name=file_name_field,
            n_peaks=basic_local["n_peaks"],
            total_bp=basic_local["total_bp"],
            min_length=basic_local["min_length"],
            max_length=basic_local["max_length"],
            mean_length=basic_local["mean_length"],
            median_length=basic_local["median_length"],
            coding_peaks=coding_peaks_local,
            coding_bp=coding_bp_local,
            coding_percent=coding_percent_local,
        )

    # Annotate per-file dataframes (if GTF provided)
    per_file_dfs = [_annotate_df(df) for df in per_file_dfs]
    # Compute per-file stats
    for p, df in zip(bed_files, per_file_dfs):
        stats_records.append(_compute_stats(df, p, p.name))

    # prepare display names for summary and processed frames
    display_names = file_names or per_file_names
    processed = process_per_file_frames(per_file_dfs)
    # If a merged bed is provided, read and prepend as the first sheet and compute stats
    sheet_names = file_names or per_file_names
    if merged_bed_path is not None and Path(merged_bed_path).exists():
        raw = read_bed_file(Path(merged_bed_path))
        # Determine if merged has 3 or 4 columns
        merged_df = raw
        if raw.shape[1] >= 4:
            # name extras with provided id column name for the 4th
            merged_df = normalize_bed_columns(raw, extra_column_names=[id_column_name] + [f"col_{i}" for i in range(5, raw.shape[1] + 1)])
        else:
            merged_df = normalize_bed_columns(raw)
        # Optionally annotate merged_df (if GTF provided)
        merged_df = _annotate_df(merged_df)
        merged_processed = process_per_file_frames([merged_df])[0]
        processed = [merged_processed] + processed
        sheet_names = [merged_sheet_name] + sheet_names
        # Compute stats for merged and prepend to stats_records
        merged_stats = _compute_stats(merged_df, Path(merged_bed_path), Path(merged_bed_path).name)
        stats_records = [merged_stats] + stats_records
        display_names = [merged_sheet_name] + display_names
    # Clean tmp files
    for t in tmp_paths:
        try:
            t.unlink(missing_ok=True)
        except Exception:
            pass
    stats_df = pd.DataFrame([s.model_dump() for s in stats_records])
    summary = compute_summary_df(stats_df, display_names=display_names,
                                 genome_total_bp=genome_total_bp, 
                                 coding_total_bp=coding_total_bp,
                                 genome_build=genome_build)
    save_excel_report(processed, sheet_names, summary, output_excel)
    typer.echo(f"Report written to: {output_excel}")


@app.command("merge-report")
def merge_report(
    bed_files: List[Path] = typer.Argument(..., exists=True, readable=True, help="BED files to analyze"),
    genome_build: str = typer.Option(..., "--genome-build", help="Genome build"),
    output: Path = typer.Option(..., "--output", "-o", help="Output Excel path"),
    refseq_gtf: Optional[Path] = typer.Option(None, "--refseq-gtf", help="RefSeq GTF for coding regions"),
    coding_bed: Optional[Path] = typer.Option(None, "--coding-bed", help="Existing or output coding BED"),
) -> None:
    """Generate an Excel report from 1..N BED files with optional coding-region analysis."""
    _generate_report(
        bed_files=bed_files,
        output_excel=output,
        refseq_gtf=refseq_gtf,
        coding_bed=coding_bed,
        genome_build=genome_build,
        extra_column_names=None,
        file_names=None,
    )


@app.command("extract-coding")
def extract_coding(
    gtf: Path = typer.Option(..., "--gtf", help="RefSeq GTF file"),
    output: Path = typer.Option(..., "--output", "-o", help="Output coding BED path"),
) -> None:
    if not check_bedtools_available():
        typer.echo("Error: bedtools not found in PATH.", err=True)
        raise typer.Exit(code=1)
    regions, total_bp = extract_coding_regions_from_gtf(gtf, output)
    typer.echo(f"Coding regions: {regions}, total bp: {total_bp:,}")


def _run_from_config(cfg: RunConfig) -> None:
    out_dir = cfg.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # Write a default report path inside the output_dir
    report_path = out_dir / "report.xlsx"

    # Optional merge step (artifact only)
    if cfg.merge.enabled:
        if not check_bedtools_available():
            typer.echo("Error: bedtools not found in PATH but required for merge.", err=True)
            raise typer.Exit(code=1)
        merged_out = cfg.merge.out_path or (out_dir / "merged_inputs.bed")
        tmp_dir = Path(cfg.tmp_dir)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        tmp_concat = tmp_dir / "_all_inputs.tmp.bed"
        # Build merge input from processed DataFrames to include ID where requested
        bed_paths = [spec.path for spec in cfg.input_config.bed_files]
        file_names = [spec.name or Path(spec.path).name for spec in cfg.input_config.bed_files]
        extra_names = [spec.extra_columns for spec in cfg.input_config.bed_files]
        per_file_dfs: List[pd.DataFrame] = []
        for idx, p in enumerate(bed_paths):
            df = load_bed_to_dataframe(p, extra_column_names=extra_names[idx])
            if cfg.input_config.add_id:
                display_name = file_names[idx]
                df = df.copy()
                df[cfg.input_config.id_column_name] = [f"{display_name}.{i+1}" for i in range(len(df))]
            per_file_dfs.append(df)
        # Write combined tmp with 3 or 4 columns depending on add_id
        with open(tmp_concat, "w") as handle:
            for idx, df in enumerate(per_file_dfs):
                if cfg.input_config.add_id and cfg.input_config.id_column_name in df.columns:
                    df[["chr", "start", "end", cfg.input_config.id_column_name]].to_csv(handle, sep="\t", header=False, index=False)
                else:
                    # Prefer an existing extra column (e.g., ID) as 4th column if present
                    preferred_extras = extra_names[idx] or []
                    selected_extra = None
                    for col_name in preferred_extras:
                        if col_name in df.columns:
                            selected_extra = col_name
                            break
                    # If no preferred name matched but there are extra columns, take the first extra
                    if selected_extra is None and len(df.columns) > 3:
                        selected_extra = df.columns[3]
                    if selected_extra is not None:
                        df[["chr", "start", "end", selected_extra]].to_csv(handle, sep="\t", header=False, index=False)
                    else:
                        df[["chr", "start", "end"]].to_csv(handle, sep="\t", header=False, index=False)
        # Build merge command
        if cfg.input_config.add_id:
            # collapse ids in 4th column
            bedtools_opts = (cfg.merge.bedtools_opts or "").strip()
            merge_cmd = [
                "bash", "-c",
                f"sort -k1,1 -k2,2n {tmp_concat} | bedtools merge -i - -c 4 -o collapse {bedtools_opts} > {merged_out}"
            ]
        else:
            sort_cmd = f"sort -k1,1 -k2,2n {tmp_concat}"
            merge_cmd = ["bash", "-c", f"{sort_cmd} | bedtools merge -i - {cfg.merge.bedtools_opts or ''} > {merged_out}"]
        subprocess = __import__("subprocess")
        subprocess.run(merge_cmd, check=True)
        tmp_concat.unlink(missing_ok=True)

    # Load genome totals and ensure coding BED exists (create if missing)
    genome_total_bp = None
    coding_total_bp = None

    with resources.files("bed_file_merger.resources").joinpath("genome_info.json").open("r") as f:
        gi = json.load(f)

    gb = gi.get(cfg.genome_build)
    if gb is None:
        typer.echo(f"Error: genome build '{cfg.genome_build}' not found in genome_info.json", err=True)
        raise typer.Exit(code=1)

    # Enforce bedtools availability
    if not check_bedtools_available():
        typer.echo("Error: bedtools not found in PATH but required for coding analysis.", err=True)
        raise typer.Exit(code=1)

    # GTF is mandatory in this pipeline
    if cfg.input_config.refseq_gtf is None or not Path(cfg.input_config.refseq_gtf).exists() or Path(cfg.input_config.refseq_gtf).stat().st_size == 0:
        typer.echo("Error: A non-empty RefSeq GTF must be provided via 'refseq_gtf' in the config.", err=True)
        raise typer.Exit(code=1)

    # Resolve coding BED path
    if cfg.input_config.coding_bed:
        coding_bed = Path(cfg.input_config.coding_bed)
    else:
        coding_bed = Path(cfg.output_dir) / f"{cfg.genome_build}_coding_regions_refseq.bed"
    coding_bed.parent.mkdir(parents=True, exist_ok=True)

    # Create coding BED if missing or empty using annotation derivation
    ann_bed = Path(cfg.output_dir) / "annotation_from_gtf.bed"
    if not ann_bed.exists() or ann_bed.stat().st_size == 0:
        build_annotation_bed_from_gtf(Path(cfg.input_config.refseq_gtf), ann_bed)
    if not coding_bed.exists() or coding_bed.stat().st_size == 0:
        derive_coding_bed_from_annotation(ann_bed, coding_bed)

    # Compute totals
    val = gb.get("Total bases")
    genome_total_bp = int(val)
    awk = sp.run(["awk", '-F\t', '{sum += $3 - $2} END {print sum}', str(coding_bed)], capture_output=True, text=True)
    try:
        coding_total_bp = int((awk.stdout or "").strip() or "0")
    except ValueError:
        coding_total_bp = 0
    if coding_total_bp <= 0:
        typer.echo("Error: coding BED creation failed or has zero total bp.", err=True)
        raise typer.Exit(code=1)

    # Generate report using config-specified display/column names
    bed_paths = [spec.path for spec in cfg.input_config.bed_files]
    file_names = [spec.name or Path(spec.path).name for spec in cfg.input_config.bed_files]
    extra_names = [spec.extra_columns for spec in cfg.input_config.bed_files]
    _generate_report(
        bed_files=bed_paths,
        output_excel=report_path,
        refseq_gtf=cfg.input_config.refseq_gtf,
        coding_bed=coding_bed,
        genome_build=cfg.genome_build,
        extra_column_names=extra_names,
        file_names=file_names,
        genome_total_bp=genome_total_bp,
        coding_total_bp=coding_total_bp,
        tmp_dir=Path(cfg.tmp_dir),
        add_id=cfg.input_config.add_id,
        id_column_name=cfg.input_config.id_column_name,
        merged_bed_path=(cfg.merge.out_path or (cfg.output_dir / "merged_inputs.bed")) if cfg.merge.enabled else None,
        merged_sheet_name="merged",
    )


@app.command("run-config")
def run_config(
    config: Path = typer.Option(..., "--config", "-c", exists=True, readable=True, help="YAML config file"),
) -> None:
    with open(config, "r") as f:
        data = yaml.safe_load(f)
    cfg = RunConfig(**data)
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    # copy the yaml into output dir for provenance
    shutil.copy2(config, cfg.output_dir / Path(config).name)
    # Optionally copy inputs for provenance
    if cfg.input_config.copy_input:
        input_copy_dir = cfg.output_dir / "input_files"
        input_copy_dir.mkdir(parents=True, exist_ok=True)

        def _deduped_destination(target_dir: Path, name: str) -> Path:
            base = Path(name).name
            candidate = target_dir / base
            if not candidate.exists():
                return candidate
            stem = Path(base).stem
            suffix = Path(base).suffix
            idx = 1
            while True:
                c = target_dir / f"{stem}_{idx}{suffix}"
                if not c.exists():
                    return c
                idx += 1

        def _is_within_output_dir(path: Path, output_dir: Path) -> bool:
            try:
                path.resolve().relative_to(output_dir.resolve())
                return True
            except Exception:
                return False

        def _copy_if_needed(src: Path, target_dir: Path) -> None:
            # Skip if already under output_dir
            if _is_within_output_dir(src, cfg.output_dir):
                return
            dest = target_dir / src.name
            if dest.exists():
                try:
                    if filecmp.cmp(src, dest, shallow=False):
                        # Identical file already present; do not copy again
                        return
                except Exception:
                    pass
                # Different content with same name: dedupe
                dest = _deduped_destination(target_dir, src.name)
            shutil.copy2(src, dest)

        # Copy bed files
        for spec in cfg.input_config.bed_files:
            _copy_if_needed(Path(spec.path), input_copy_dir)
        # Copy GTF if provided
        if cfg.input_config.refseq_gtf is not None:
            _copy_if_needed(Path(cfg.input_config.refseq_gtf), input_copy_dir)
    _run_from_config(cfg)


def main() -> None:
    app()


if __name__ == "__main__":
    main()


