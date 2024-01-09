"""Definition of the PRP command-line interface."""
import json
import logging
from typing import List

import click
import pandas as pd
from pydantic import TypeAdapter, ValidationError

from .models.metadata import SoupType, SoupVersion
from .models.phenotype import ElementType
from .models.qc import QcMethodIndex, QcSoftware
from .models.sample import MethodIndex, PipelineResult
from .models.typing import TypingMethod
from .parse import (
    parse_amrfinder_amr_pred,
    parse_amrfinder_vir_pred,
    parse_cgmlst_results,
    parse_kraken_result,
    parse_mlst_results,
    parse_mykrobe_amr_pred,
    parse_mykrobe_lineage_results,
    parse_postalignqc_results,
    parse_quast_results,
    parse_resfinder_amr_pred,
    parse_tbprofiler_amr_pred,
    parse_tbprofiler_lineage_results,
    parse_virulencefinder_stx_typing,
    parse_virulencefinder_vir_pred,
)
from .parse.metadata import get_database_info, parse_run_info

logging.basicConfig(
    level=logging.INFO, format="[%(asctime)s] %(levelname)s in %(module)s: %(message)s"
)
LOG = logging.getLogger(__name__)

OUTPUT_SCHEMA_VERSION = 1


@click.group()
def cli():
    """Base CLI entrypoint."""


@cli.command()
@click.option("-i", "--sample-id", required=True, help="Sample identifier")
@click.option(
    "-u",
    "--run-metadata",
    type=click.File(),
    required=True,
    help="Analysis metadata from the pipeline in json format",
)
@click.option("-q", "--quast", type=click.File(), help="Quast quality control metrics")
@click.option(
    "-d",
    "--process-metadata",
    type=click.File(),
    multiple=True,
    help="Nextflow processes metadata from the pipeline in json format",
)
@click.option(
    "-k", "--kraken", type=click.File(), help="Kraken species annotation results"
)
@click.option(
    "-a",
    "--amrfinder",
    type=click.Path(),
    help="amrfinderplus anti-microbial resistance results",
)
@click.option("-m", "--mlst", type=click.File(), help="MLST prediction results")
@click.option("-c", "--cgmlst", type=click.File(), help="cgMLST prediction results")
@click.option(
    "-v",
    "--virulencefinder",
    type=click.Path(),
    help="Virulence factor prediction results",
)
@click.option(
    "-r",
    "--resfinder",
    type=click.File(),
    help="resfinder resistance prediction results",
)
@click.option("-p", "--quality", type=click.File(), help="postalignqc qc results")
@click.option("-k", "--mykrobe", type=click.File(), help="mykrobe results")
@click.option("-t", "--tbprofiler", type=click.File(), help="tbprofiler results")
@click.option("--correct_alleles", is_flag=True, help="Correct alleles")
@click.option(
    "-o", "--output", required=True, type=click.File("w"), help="output filepath"
)
def create_bonsai_input(
    sample_id,
    run_metadata,
    quast,
    process_metadata,
    kraken,
    mlst,
    cgmlst,
    virulencefinder,
    amrfinder,
    resfinder,
    quality,
    mykrobe,
    tbprofiler,
    correct_alleles,
    output,
):  # pylint: disable=too-many-arguments
    """Combine pipeline results into a standardized json output file."""
    LOG.info("Start generating pipeline result json")
    results = {
        "run_metadata": {
            "run": parse_run_info(run_metadata),
            "databases": get_database_info(process_metadata),
        },
        "qc": [],
        "typing_result": [],
        "element_type_result": [],
    }
    # qc
    if quast:
        LOG.info("Parse quast results")
        res: QcMethodIndex = parse_quast_results(quast)
        results["qc"].append(res)
    if quality:
        LOG.info("Parse quality results")
        res: QcMethodIndex = parse_postalignqc_results(quality)
        results["qc"].append(res)

    # typing
    if mlst:
        LOG.info("Parse mlst results")
        res: MethodIndex = parse_mlst_results(mlst)
        results["typing_result"].append(res)
    if cgmlst:
        LOG.info("Parse cgmlst results")
        res: MethodIndex = parse_cgmlst_results(cgmlst, correct_alleles=correct_alleles)
        results["typing_result"].append(res)

    # resfinder of different types
    if resfinder:
        LOG.info("Parse resistance results")
        pred_res = json.load(resfinder)
        methods = [
            ElementType.AMR,
            ElementType.STRESS,
        ]
        for method in methods:
            res: MethodIndex = parse_resfinder_amr_pred(pred_res, method)
            # exclude empty results from output
            if len(res.result.genes) > 0 and len(res.result.variants) > 0:
                results["element_type_result"].append(res)

    # amrfinder
    if amrfinder:
        LOG.info("Parse amr results")
        methods = [
            ElementType.AMR,
            ElementType.STRESS,
        ]
        for method in methods:
            res: MethodIndex = parse_amrfinder_amr_pred(amrfinder, method)
            results["element_type_result"].append(res)
        vir: MethodIndex = parse_amrfinder_vir_pred(amrfinder)
        results["element_type_result"].append(vir)

    # get virulence factors in sample
    if virulencefinder:
        LOG.info("Parse virulencefinder results")
        # virulence genes
        vir: MethodIndex | None = parse_virulencefinder_vir_pred(virulencefinder)
        if vir is not None:
            results["element_type_result"].append(vir)

        # stx typing
        res: MethodIndex | None = parse_virulencefinder_stx_typing(virulencefinder)
        if res is not None:
            results["typing_result"].append(res)

    # species id
    if kraken:
        LOG.info("Parse kraken results")
        results["species_prediction"] = parse_kraken_result(kraken)
    else:
        results["species_prediction"] = []

    # mycobacterium tuberculosis
    # mykrobe
    if mykrobe:
        LOG.info("Parse mykrobe results")
        pred_res = pd.read_csv(mykrobe, quotechar='"')
        pred_res.columns.values[3] = "variants"
        pred_res.columns.values[4] = "genes"
        pred_res = pred_res.to_dict(orient="records")

        # verify that sample id is in prediction result
        if not sample_id in pred_res[0]["sample"]:
            LOG.warning(
                "Sample id %s is not in Mykrobe result, possible sample mixup",
                sample_id,
            )
            raise click.Abort()

        # add mykrobe db version
        results["run_metadata"]["databases"].append(
            SoupVersion(
                name="mykrobe-predictor",
                version=pred_res[0]["mykrobe_version"],
                type=SoupType.DB,
            )
        )
        # parse mykrobe result
        amr_res = parse_mykrobe_amr_pred(pred_res, ElementType.AMR)
        if amr_res is not None:
            results["element_type_result"].append(amr_res)

        lin_res: MethodIndex | None = parse_mykrobe_lineage_results(
            pred_res, TypingMethod.LINEAGE
        )
        if lin_res is not None:
            results["typing_result"].append(lin_res)

    # tbprofiler
    if tbprofiler:
        LOG.info("Parse tbprofiler results")
        pred_res = json.load(tbprofiler)
        db_info: List[SoupVersion] = []
        db_info = [
            SoupVersion(
                name=pred_res["db_version"]["name"],
                version=pred_res["db_version"]["commit"],
                type=SoupType.DB,
            )
        ]
        results["run_metadata"]["databases"].extend(db_info)
        lin_res: MethodIndex = parse_tbprofiler_lineage_results(
            pred_res, TypingMethod.LINEAGE
        )
        results["typing_result"].append(lin_res)
        amr_res: MethodIndex = parse_tbprofiler_amr_pred(pred_res, ElementType.AMR)
        results["element_type_result"].append(amr_res)

    try:
        output_data = PipelineResult(
            sample_id=sample_id, schema_version=OUTPUT_SCHEMA_VERSION, **results
        )
    except ValidationError as err:
        click.secho("Input failed Validation", fg="red")
        click.secho(err)
        raise click.Abort
    LOG.info("Storing results to: %s", output.name)
    output.write(output_data.model_dump_json(indent=2))
    click.secho("Finished generating pipeline output", fg="green")


@cli.command()
def print_schema():
    """Print Pipeline result output format schema."""
    click.secho(PipelineResult.schema_json(indent=2))


@cli.command()
@click.option("-o", "--output", required=True, type=click.File("r"))
def validate(output):
    """Validate output format of result json file."""
    js = json.load(output)
    try:
        PipelineResult(**js)
    except ValidationError as err:
        click.secho("Invalid file format X", fg="red")
        click.secho(err)
    else:
        click.secho(f'The file "{output.name}" is valid', fg="green")


@cli.command()
@click.option("-q", "--quast", type=click.File(), help="Quast quality control metrics")
@click.option("-p", "--quality", type=click.File(), help="postalignqc qc results")
@click.option("-c", "--cgmlst", type=click.File(), help="cgMLST prediction results")
@click.option("--correct_alleles", is_flag=True, help="Correct alleles")
@click.option(
    "-o", "--output", required=True, type=click.File("w"), help="output filepath"
)
def create_cdm_input(quast, quality, cgmlst, correct_alleles, output) -> None:
    """Format QC metrics into CDM compatible input file."""
    results = []
    if quality:
        LOG.info("Parse quality results")
        res: QcMethodIndex = parse_postalignqc_results(quality)
        results.append(res)

    if quast:
        LOG.info("Parse quast results")
        res: QcMethodIndex = parse_quast_results(quast)
        results.append(res)

    if cgmlst:
        LOG.info("Parse cgmlst results")
        res: MethodIndex = parse_cgmlst_results(cgmlst, correct_alleles=correct_alleles)
        n_missing_loci = QcMethodIndex(
            software=QcSoftware.CHEWBBACA, result={"n_missing": res.result.n_missing}
        )
        results.append(n_missing_loci)
    # cast output as pydantic type for easy serialization
    qc_data = TypeAdapter(List[QcMethodIndex])

    LOG.info("Storing results to: %s", output.name)
    output.write(qc_data.dump_json(results, indent=3).decode("utf-8"))
    click.secho("Finished generating QC output", fg="green")
