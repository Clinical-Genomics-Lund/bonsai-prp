"""Shared utility functions."""
from typing import Tuple

from ...models.phenotype import (
    ElementType,
    ElementTypeResult,
    PhenotypeInfo,
    VariantType,
)


def _default_amr_phenotype() -> PhenotypeInfo:
    return PhenotypeInfo(
        type=ElementType.AMR,
        group=ElementType.AMR,
        name=ElementType.AMR,
    )


def is_prediction_result_empty(result: ElementTypeResult) -> bool:
    """Check if prediction result is emtpy.

    :param result: Prediction result
    :type result: ElementTypeResult
    :return: Retrun True if no resistance was predicted.
    :rtype: bool
    """
    n_entries = len(result.genes) + len(result.variants)
    return n_entries == 0


def get_nt_change(ref_codon: str, alt_codon: str) -> Tuple[str, str]:
    """Get nucleotide change from codons

    Ref: TCG, Alt: TTG => Tuple[C, T]

    :param ref_codon: Reference codeon
    :type ref_codon: str
    :param str: Alternatve codon
    :type str: str
    :return: Returns nucleotide changed from the reference.
    :rtype: Tuple[str, str]
    """
    ref_nt = ""
    alt_nt = ""
    for ref, alt in zip(ref_codon, alt_codon):
        if not ref == alt:
            ref_nt += ref
            alt_nt += alt
    return ref_nt.upper(), alt_nt.upper()


def format_nt_change(
    ref: str,
    alt: str,
    var_type: VariantType,
    start_pos: int,
    end_pos: int = None,
) -> str:
    """Format nucleotide change

    :param ref: Reference sequence
    :type ref: str
    :param alt: Alternate sequence
    :type alt: str
    :param pos: Position
    :type pos: int
    :param var_type: Type of change
    :type var_type: VariantType
    :return: Formatted nucleotide
    :rtype: str
    """
    fmt_change = ""
    match var_type:
        case VariantType.SUBSTITUTION:
            f"g.{start_pos}{ref}>{alt}"
        case VariantType.DELETION:
            f"g.{start_pos}_{end_pos}del"
        case VariantType.INSERTION:
            f"g.{start_pos}_{end_pos}ins{alt}"
    return fmt_change
