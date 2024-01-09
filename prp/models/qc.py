"""QC data models."""
from enum import Enum
from typing import Dict

from pydantic import BaseModel, Field

from .base import RWModel
from .typing import TypingSoftware


class QcSoftware(Enum):
    """Valid tools."""

    QUAST = "quast"
    FASTQC = "fastqc"
    POSTALIGNQC = "postalignqc"
    CHEWBBACA = TypingSoftware.CHEWBBACA.value


class QuastQcResult(BaseModel):
    """Assembly QC metrics."""

    total_length: int
    reference_length: int
    largest_contig: int
    n_contigs: int
    n50: int
    assembly_gc: float
    reference_gc: float
    duplication_ratio: float


class PostAlignQcResult(BaseModel):
    """Alignment QC metrics."""

    ins_size: int
    ins_size_dev: int
    mean_cov: int
    pct_above_x: Dict[str, float]
    mapped_reads: int
    tot_reads: int
    iqr_median: float


class GenomeCompleteness(BaseModel):
    """Alignment QC metrics."""

    n_missing: int = Field(..., description="Number of missing cgMLST alleles")


class QcMethodIndex(RWModel):
    """QC results container.

    Based on Mongo db Attribute pattern.
    Reference: https://www.mongodb.com/developer/products/mongodb/attribute-pattern/
    """

    software: QcSoftware
    version: str | None = None
    result: QuastQcResult | PostAlignQcResult | GenomeCompleteness
