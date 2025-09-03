from __future__ import annotations

"""Pydantic models for BED analysis configuration and results."""

from pathlib import Path
from typing import List, Optional

from pydantic import BaseModel, Field


class AnalysisConfig(BaseModel):
    bed_files: List[Path] = Field(..., description="BED files to analyze")
    output_excel: Path = Field(..., description="Path to output Excel report")
    refseq_gtf: Optional[Path] = Field(None, description="RefSeq GTF for coding regions")
    coding_bed: Optional[Path] = Field(None, description="Precomputed coding regions BED")


class BedFileStats(BaseModel):
    file_name: str
    n_peaks: int
    total_bp: int
    min_length: int
    max_length: int
    mean_length: float
    median_length: float
    coding_peaks: Optional[int] = None
    coding_bp: Optional[int] = None
    coding_percent: Optional[float] = None


class MergeOptions(BaseModel):
    enabled: bool = Field(default=False, description="Whether to run bedtools merge on all inputs")
    out_path: Optional[Path] = Field(None, description="Path to write merged BED if enabled")
    bedtools_opts: Optional[str] = Field(None, description="Extra bedtools merge options")


class RunConfig(BaseModel):
    bed_files: List[Path]
    output_dir: Path
    refseq_gtf: Optional[Path] = None
    coding_bed: Optional[Path] = None
    merge: MergeOptions = Field(default_factory=MergeOptions)
    file_names: Optional[List[str]] = Field(None, description="Optional display names for each BED file")
    extra_column_names: Optional[List[Optional[List[str]]]] = Field(
        None,
        description="Optional per-file list of extra column names (4th+ columns)"
    )
    genome_build: Optional[str] = Field(None, description="Genome build key, e.g., hg38 or mm10")


