"""Datamodels used for prediction results."""
from enum import Enum
from typing import Any, Dict, List, Optional, Union

from pydantic import BaseModel, Field

from .base import RWModel


class SequenceStand(Enum):
    """Definition of DNA strand."""

    FORWARD = "+"
    REVERSE = "-"


class PredictionSoftware(Enum):
    """Container for prediciton software names."""

    AMRFINDER = "amrfinder"
    RESFINDER = "resfinder"
    VIRFINDER = "virulencefinder"
    MYKROBE = "mykrobe"
    TBPROFILER = "tbprofiler"


class VariantType(Enum):
    """Types of variants."""

    SUBSTITUTION = "substitution"
    INSERTION = "insertion"
    DELETION = "deletion"


class ElementType(Enum):
    """Categories of resistance and virulence genes."""

    AMR = "AMR"
    STRESS = "STRESS"
    VIR = "VIRULENCE"


class ElementStressSubtype(Enum):
    """Categories of resistance and virulence genes."""

    ACID = "ACID"
    BIOCIDE = "BIOCIDE"
    METAL = "METAL"
    HEAT = "HEAT"


class ElementAmrSubtype(Enum):
    """Categories of resistance and virulence genes."""

    AMR = "AMR"


class ElementVirulenceSubtype(Enum):
    """Categories of resistance and virulence genes."""

    VIR = "VIRULENCE"
    ANTIGEN = "ANTIGEN"
    TOXIN = "TOXIN"


class PhenotypeInfo(RWModel):
    """Phenotype information."""

    name: str
    group: str | None = Field(None, description="Name of the group a trait belongs to.")
    type: ElementType = Field(
        ..., description="Trait category, for example AMR, STRESS etc."
    )
    reference: List[str] = Field([], description="References supporting trait")
    note: str | None = Field(None, description="Note, can be used for confidence score")


class DatabaseReference(RWModel):
    """Refernece to a database."""

    ref_database: Optional[str] = None
    ref_id: Optional[str] = None


class GeneBase(BaseModel):
    """Container for gene information"""

    # basic info
    gene_symbol: Optional[str] = None
    accession: Optional[str] = None
    sequence_name: Optional[str] = Field(
        default=None, description="Reference sequence name"
    )
    element_type: ElementType = Field(
        description="The predominant function fo the gene."
    )
    element_subtype: Union[
        ElementStressSubtype, ElementAmrSubtype, ElementVirulenceSubtype
    ] = Field(description="Further functional categorization of the genes.")
    # position
    ref_start_pos: Optional[int] = Field(
        None, description="Alignment start in reference"
    )
    ref_end_pos: Optional[int] = Field(None, description="Alignment end in reference")
    ref_gene_length: Optional[int] = Field(
        default=None,
        alias="target_length",
        description="The length of the reference protein or gene.",
    )

    # prediction
    method: Optional[str] = Field(None, description="Method used to predict gene")
    identity: Optional[float] = Field(
        None, description="Identity to reference sequence"
    )
    coverage: Optional[float] = Field(
        None, description="Ratio reference sequence covered"
    )


class AmrFinderGene(GeneBase):
    """Container for Resfinder gene prediction information"""

    contig_id: str
    query_start_pos: int = Field(
        default=None, description="Start position on the assembly"
    )
    query_end_pos: int = Field(default=None, description="End position on the assembly")
    strand: SequenceStand


class AmrFinderResistanceGene(AmrFinderGene):
    """AMRfinder resistance gene information."""

    phenotypes: List[PhenotypeInfo] = []


class ResistanceGene(GeneBase):
    """Container for resistance gene information"""

    phenotypes: List[PhenotypeInfo] = []


class VirulenceGene(GeneBase, DatabaseReference):
    """Container for virulence gene information"""

    depth: Optional[float] = Field(
        None, description="Ammount of sequence data supporting the gene."
    )


class ResfinderGene(ResistanceGene):
    """Container for Resfinder gene prediction information"""

    depth: Optional[float] = Field(
        None, description="Ammount of sequence data supporting the gene."
    )


class VariantBase(RWModel):
    """Container for mutation information"""

    # classification
    variant_type: VariantType
    phenotypes: List[PhenotypeInfo] = []

    # variant location
    gene_symbol: str
    accession: Optional[str] = None
    position: int
    ref_nt: str
    alt_nt: str
    ref_aa: Optional[str] = None
    alt_aa: Optional[str] = None

    # prediction info
    depth: Optional[float] = Field(None, description="Total depth, ref + alt.")
    frequency: Optional[float] = Field(None, description="Alt allele frequency.")
    confidence: Optional[int] = Field(None, description="Genotype confidence.")
    method: str = Field(..., description="Prediction method used to call variant")
    passed_qc: bool = Field(
        ..., description="Describe if variant has passed the tool qc check"
    )


class ResfinderVariant(VariantBase):
    """Container for ResFinder variant information"""


class MykrobeVariant(VariantBase):
    """Container for Mykrobe variant information"""


class TbProfilerVariant(VariantBase):
    """Container for TbProfiler variant information"""

    variant_effect: str
    hgvs_nt_change: Optional[str] = Field(..., description="DNA change in HGVS format")
    hgvs_aa_change: Optional[str] = Field(
        ..., description="Protein change in HGVS format"
    )


class ElementTypeResult(BaseModel):
    """Phenotype result data model.

    A phenotype result is a generic data structure that stores predicted genes,
    mutations and phenotyp changes.
    """

    phenotypes: Dict[str, List[str]]
    genes: List[Union[AmrFinderResistanceGene, AmrFinderGene, ResfinderGene, VirulenceGene]]
    variants: List[Union[TbProfilerVariant, MykrobeVariant, ResfinderVariant]]