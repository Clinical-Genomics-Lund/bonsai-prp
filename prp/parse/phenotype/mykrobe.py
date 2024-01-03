"""Parse Mykrobe results."""
import logging
import re
from typing import Any, Dict, Tuple

from ...models.phenotype import ElementAmrSubtype, ElementType, ElementTypeResult
from ...models.phenotype import PredictionSoftware as Software
from ...models.phenotype import ResistanceGene, ResistanceVariant, VariantType
from ...models.sample import MethodIndex
from .utils import is_prediction_result_empty

LOG = logging.getLogger(__name__)


def _get_mykrobe_amr_sr_profie(mykrobe_result):
    """Get mykrobe susceptibility/resistance profile."""
    susceptible = set()
    resistant = set()

    if not mykrobe_result:
        return {}

    for element_type in mykrobe_result:
        if element_type["susceptibility"].upper() == "R":
            resistant.add(element_type["drug"])
        else:
            susceptible.add(element_type["drug"])
    return {"susceptible": list(susceptible), "resistant": list(resistant)}


def get_mutation_type(var_nom: str) -> Tuple[VariantType, str, str, int]:
    """Extract mutation type from Mykrobe mutation description.

    GCG7569GTG -> mutation type, ref_nt, alt_nt, pos

    :param var_nom: Mykrobe mutation description
    :type var_nom: str
    :return: Return variant type, ref_codon, alt_codont and position
    :rtype: Tuple[VariantType, str, str, int]
    """
    mut_type = None
    ref_codon = None
    alt_codon = None
    position = None
    try:
        ref_idx = re.search(r"\d", var_nom, 1).start()
        alt_idx = re.search(r"\d(?=[^\d]*$)", var_nom).start() + 1
    except AttributeError:
        return mut_type, ref_codon, alt_codon, position

    ref_codon = var_nom[:ref_idx]
    alt_codon = var_nom[alt_idx:]
    position = int(var_nom[ref_idx:alt_idx])
    if len(ref_codon) > len(alt_codon):
        mut_type = VariantType.DELETION
    elif len(ref_codon) < len(alt_codon):
        mut_type = VariantType.INSERTION
    else:
        mut_type = VariantType.SUBSTITUTION
    return mut_type, ref_codon, alt_codon, position


def _parse_mykrobe_amr_variants(mykrobe_result) -> Tuple[ResistanceVariant, ...]:
    """Get resistance genes from mykrobe result."""
    results = []

    for element_type in mykrobe_result:
        # skip non-resistance yeilding
        if not element_type["susceptibility"].upper() == "R":
            continue

        if element_type["variants"] is None:
            continue

        try:
            depth = float(element_type["genes"].split(':')[-1])
        except AttributeError:
            depth = None

        var_info = element_type["variants"].split("-")[1].split(":")[0]
        _, ref_nt, alt_nt, position = get_mutation_type(var_info)
        var_nom = element_type["variants"].split("-")[0].split("_")[1]
        var_type, *_ = get_mutation_type(var_nom)
        variant = ResistanceVariant(
            variant_type=var_type,
            gene_symbol=element_type["variants"].split("_")[0],
            position=position,
            ref_nt=ref_nt,
            alt_nt=alt_nt,
            depth=depth,
            change=var_nom,
            drugs=[element_type["drug"].lower()],
        )
        results.append(variant)
    return results


def parse_mykrobe_amr_pred(
    prediction: Dict[str, Any], resistance_category
) -> ElementTypeResult | None:
    """Parse mykrobe resistance prediction results."""
    LOG.info("Parsing mykrobe prediction")
    resistance = ElementTypeResult(
        phenotypes=_get_mykrobe_amr_sr_profie(prediction),
        genes=[],
        mutations=_parse_mykrobe_amr_variants(prediction),
    )

    # verify prediction result
    if is_prediction_result_empty(resistance):
        result = None
    else:
        result = MethodIndex(
            type=resistance_category, software=Software.MYKROBE, result=resistance
        )
    return result
