"""Parse TBprofiler result."""
import logging
from typing import Any, Dict, List, Tuple

from ...models.metadata import SoupVersions
from ...models.phenotype import ElementType, ElementTypeResult, PhenotypeInfo
from ...models.phenotype import PredictionSoftware as Software
from ...models.phenotype import TbProfilerVariant, VariantType
from ...models.sample import MethodIndex

LOG = logging.getLogger(__name__)


def _get_tbprofiler_amr_sr_profie(tbprofiler_result):
    """Get tbprofiler susceptibility/resistance profile."""
    susceptible = set()
    resistant = set()
    drugs = [
        "ofloxacin",
        "moxifloxacin",
        "isoniazid",
        "delamanid",
        "kanamycin",
        "amikacin",
        "ethambutol",
        "ethionamide",
        "streptomycin",
        "ciprofloxacin",
        "levofloxacin",
        "pyrazinamide",
        "linezolid",
        "rifampicin",
        "capreomycin",
    ]

    if not tbprofiler_result:
        return {}

    for hit in tbprofiler_result["dr_variants"]:
        for drug in hit["gene_associated_drugs"]:
            resistant.add(drug)
    susceptible = [drug for drug in drugs if drug not in resistant]
    return {"susceptible": list(susceptible), "resistant": list(resistant)}


def _parse_tbprofiler_amr_variants(predictions) -> Tuple[TbProfilerVariant, ...]:
    """Get resistance genes from tbprofiler result."""
    variant_caller = None
    for prog in predictions["pipeline"]:
        if prog["Analysis"].lower() == "variant calling":
            variant_caller = prog["Program"]
    results = []

    # tbprofiler report three categories of variants
    # - dr_variants: known resistance variants
    # - qc_fail_variants: known resistance variants failing qc
    # - other_variants: variants not in the database but in genes associated with resistance
    for result_type in ["dr_variants", "other_variants", "qc_fail_variants"]:
        # associated with passed/ failed qc
        if result_type == "qc_fail_variants":
            passed_qc = False
        else:
            passed_qc = True

        # parse variants
        for hit in predictions.get(result_type, []):
            ref_nt = hit["ref"]
            alt_nt = hit["alt"]
            if len(ref_nt) == len(alt_nt):
                var_type = VariantType.SUBSTITUTION
            elif len(ref_nt) > len(alt_nt):
                var_type = VariantType.DELETION
            else:
                var_type = VariantType.INSERTION

            variant = TbProfilerVariant(
                # classificatoin
                variant_type=var_type,
                phenotypes=parse_drug_resistance_info(hit.get("drugs", [])),
                # location
                gene_symbol=hit["gene"],
                accession=hit["feature_id"],
                position=int(hit["genome_pos"]),
                ref_nt=ref_nt,
                alt_nt=alt_nt,
                # consequense
                variant_effect=hit["type"],
                hgvs_nt_change=hit["nucleotide_change"],
                hgvs_aa_change=hit["protein_change"],
                # prediction info
                depth=hit["depth"],
                frequency=float(hit["freq"]),
                method=variant_caller,
                passed_qc=passed_qc,
            )
            results.append(variant)
    # sort variants
    variants = sorted(results, key=lambda entry: (entry.gene_symbol, entry.position))
    return variants


def parse_drug_resistance_info(drugs: List[Dict[str, str]]) -> List[PhenotypeInfo]:
    """Parse drug info into the standard format.

    :param drugs: TbProfiler drug info
    :type drugs: List[Dict[str, str]]
    :return: Formatted phenotype info
    :rtype: List[PhenotypeInfo]
    """
    phenotypes = []
    for drug in drugs:
        # assign element type
        if drug["type"] == "drug" and drug["confers"] == "resistance":
            drug_type = ElementType.AMR
        else:
            drug_type = ElementType.AMR
            LOG.warning(
                "Unknown TbProfiler drug; drug: %s, confers: %s; default to %s",
                drug["type"],
                drug["confers"],
                drug_type,
            )
        reference = drug.get("literature")
        phenotypes.append(
            PhenotypeInfo(
                name=drug["drug"],
                type=drug_type,
                reference=[] if reference is None else [reference],
                note=drug.get("who confidence"),
            )
        )
    return phenotypes


def parse_tbprofiler_amr_pred(
    prediction: Dict[str, Any], resistance_category
) -> Tuple[SoupVersions, ElementTypeResult]:
    """Parse tbprofiler resistance prediction results."""
    LOG.info("Parsing tbprofiler prediction")
    resistance = ElementTypeResult(
        phenotypes=_get_tbprofiler_amr_sr_profie(prediction),
        genes=[],
        variants=_parse_tbprofiler_amr_variants(prediction),
    )
    return MethodIndex(
        type=resistance_category, software=Software.TBPROFILER, result=resistance
    )
