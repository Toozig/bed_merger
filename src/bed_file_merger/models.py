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


class BedFileSpec(BaseModel):
    """Specification of a single BED input, including metadata for reporting and parsing."""

    path: Path
    name: Optional[str] = Field(default=None, description="Display name for this BED in the report")
    extra_columns: Optional[List[str]] = Field(
        default=None,
        description="Names for 4th+ columns of this BED; omit or set null to auto-name",
    )


class InputConfig(BaseModel):
    bed_files: List[BedFileSpec]
    refseq_gtf: Optional[Path] = None
    coding_bed: Optional[Path] = None
    copy_input: bool = Field(default=False, description="If true, copy inputs into output_dir/input_files")
    add_id: bool = Field(default=False, description="If true, add an id column to each peak")
    id_column_name: str = Field(default="id", description="Column name for generated peak IDs")


class MergeOptions(BaseModel):
    enabled: bool = Field(default=False, description="Whether to run bedtools merge on all inputs")
    out_path: Optional[Path] = Field(None, description="Path to write merged BED if enabled")
    bedtools_opts: Optional[str] = Field(None, description="Extra bedtools merge options")


class RunConfig(BaseModel):
    input_config: InputConfig
    output_dir: Path
    merge: MergeOptions = Field(default_factory=MergeOptions)
    genome_build: Optional[str] = Field(None, description="Genome build key, e.g., hg38 or mm10")
    tmp_dir: Path = Field(default=Path("/tmp"), description="Directory for temporary files; defaults to /tmp")


