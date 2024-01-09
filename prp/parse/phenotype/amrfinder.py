"""Parse AMRfinder plus result."""
import logging
from typing import Tuple

import pandas as pd

from ...models.phenotype import (
    AmrFinderGene,
    AmrFinderResistanceGene,
    ElementType,
    ElementTypeResult,
    PhenotypeInfo,
)
from ...models.phenotype import PredictionSoftware as Software
from ...models.sample import MethodIndex

LOG = logging.getLogger(__name__)


def _parse_amrfinder_amr_results(
    predictions: dict,
) -> Tuple[AmrFinderResistanceGene, ...]:
    """Parse amrfinder prediction results from amrfinderplus."""
    genes = []
    for prediction in predictions:
        element_type = ElementType(prediction["element_type"])
        res_class = prediction["Class"]
        res_sub_class = prediction["Subclass"]

        # classification to phenotype object
        phenotypes = []
        if res_class is None:
            phenotypes.append(
                PhenotypeInfo(
                    type=element_type,
                    group=element_type,
                    name=element_type,
                )
            )
        elif isinstance(res_sub_class, str):
            phenotypes.extend(
                [
                    PhenotypeInfo(
                        type=element_type,
                        group=res_class.lower(),
                        name=annot.lower(),
                    )
                    for annot in res_sub_class.split("/")
                ]
            )
        # store resistance gene
        gene = AmrFinderResistanceGene(
            # info
            gene_symbol=prediction["gene_symbol"],
            accession=prediction["close_seq_accn"],
            sequence_name=prediction["sequence_name"],
            # gene classification
            element_type=element_type,
            element_subtype=prediction["element_subtype"],
            phenotypes=phenotypes,
            # position
            contig_id=prediction["contig_id"],
            query_start_pos=prediction["Start"],
            query_end_pos=prediction["Stop"],
            strand=prediction["Strand"],
            ref_gene_length=prediction["ref_seq_len"],
            alignment_length=prediction["align_len"],
            # prediction
            method=prediction["Method"],
            identity=prediction["ref_seq_identity"],
            coverage=prediction["ref_seq_cov"],
        )
        genes.append(gene)

    # concat resistance profile
    sr_profile = {
        "susceptible": [],
        "resistant": list({pheno.name for gene in genes for pheno in gene.phenotypes}),
    }
    # sort genes
    genes = sorted(genes, key=lambda entry: (entry.gene_symbol, entry.coverage))
    return ElementTypeResult(phenotypes=sr_profile, genes=genes, variants=[])


def parse_amrfinder_amr_pred(file: str, element_type: ElementType) -> ElementTypeResult:
    """Parse amrfinder resistance prediction results."""
    LOG.info("Parsing amrfinder amr prediction")
    with open(file, "rb") as tsvfile:
        hits = pd.read_csv(tsvfile, delimiter="\t")
        hits = hits.rename(
            columns={
                "Contig id": "contig_id",
                "Gene symbol": "gene_symbol",
                "Sequence name": "sequence_name",
                "Element type": "element_type",
                "Element subtype": "element_subtype",
                "Target length": "target_length",
                "Reference sequence length": "ref_seq_len",
                "% Coverage of reference sequence": "ref_seq_cov",
                "% Identity to reference sequence": "ref_seq_identity",
                "Alignment length": "align_len",
                "Accession of closest sequence": "close_seq_accn",
                "Name of closest sequence": "close_seq_name",
            }
        )
        hits = hits.drop(columns=["Protein identifier", "HMM id", "HMM description"])
        hits = hits.where(pd.notnull(hits), None)
        # group predictions based on their element type
        predictions = hits.loc[
            lambda row: row.element_type == element_type.value
        ].to_dict(orient="records")
        results: ElementTypeResult = _parse_amrfinder_amr_results(predictions)
    return MethodIndex(type=element_type, result=results, software=Software.AMRFINDER)


def _parse_amrfinder_vir_results(predictions: dict) -> ElementTypeResult:
    """Parse amrfinder prediction results from amrfinderplus."""
    genes = []
    for prediction in predictions:
        gene = AmrFinderGene(
            # info
            gene_symbol=prediction["gene_symbol"],
            accession=prediction["close_seq_accn"],
            sequence_name=prediction["sequence_name"],
            # gene classification
            element_type=prediction["element_type"],
            element_subtype=prediction["element_subtype"],
            # position
            contig_id=prediction["contig_id"],
            query_start_pos=prediction["Start"],
            query_end_pos=prediction["Stop"],
            strand=prediction["Strand"],
            ref_gene_length=prediction["ref_seq_len"],
            alignment_length=prediction["align_len"],
            # prediction
            method=prediction["Method"],
            identity=prediction["ref_seq_identity"],
            coverage=prediction["ref_seq_cov"],
        )
        genes.append(gene)
    # sort genes
    genes = sorted(genes, key=lambda entry: (entry.gene_symbol, entry.coverage))
    return ElementTypeResult(phenotypes={}, genes=genes, variants=[])


def parse_amrfinder_vir_pred(file: str):
    """Parse amrfinder virulence prediction results."""
    LOG.info("Parsing amrfinder virulence prediction")
    with open(file, "rb") as tsvfile:
        hits = pd.read_csv(tsvfile, delimiter="\t")
        hits = hits.rename(
            columns={
                "Contig id": "contig_id",
                "Gene symbol": "gene_symbol",
                "Sequence name": "sequence_name",
                "Element type": "element_type",
                "Element subtype": "element_subtype",
                "Target length": "target_length",
                "Reference sequence length": "ref_seq_len",
                "% Coverage of reference sequence": "ref_seq_cov",
                "% Identity to reference sequence": "ref_seq_identity",
                "Alignment length": "align_len",
                "Accession of closest sequence": "close_seq_accn",
                "Name of closest sequence": "close_seq_name",
            }
        )
        hits = hits.drop(columns=["Protein identifier", "HMM id", "HMM description"])
        hits = hits.where(pd.notnull(hits), None)
        predictions = hits[hits["element_type"] == "VIRULENCE"].to_dict(
            orient="records"
        )
        results: ElementTypeResult = _parse_amrfinder_vir_results(predictions)
    return MethodIndex(
        type=ElementType.VIR, software=Software.AMRFINDER, result=results
    )
