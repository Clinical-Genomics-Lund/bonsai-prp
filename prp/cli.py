"""Definition of the PRP command-line interface."""
import json
import logging
from pathlib import Path
from typing import List

import click
import pandas as pd
import pysam
from cyvcf2 import VCF, Writer
from pydantic import TypeAdapter, ValidationError

from .models.metadata import SoupType, SoupVersion
from .models.phenotype import ElementType
from .models.qc import QcMethodIndex, QcSoftware
from .models.sample import MethodIndex, PipelineResult, ReferenceGenome
from .models.typing import TypingMethod
from .parse import (
    load_variants,
    parse_alignment_results,
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
    parse_serotypefinder_oh_typing,
    parse_tbprofiler_amr_pred,
    parse_tbprofiler_lineage_results,
    parse_virulencefinder_stx_typing,
    parse_virulencefinder_vir_pred,
)
from .parse.mapping import get_reference_seq_accnr
from .parse.metadata import get_database_info, get_gb_genome_version, parse_run_info
from .parse.utils import get_db_version
from .parse.variant import annotate_delly_variants

logging.basicConfig(
    level=logging.INFO, format="[%(asctime)s] %(levelname)s in %(module)s: %(message)s"
)
LOG = logging.getLogger(__name__)

OUTPUT_SCHEMA_VERSION = 1


@click.group()
def cli():
    """Jasen pipeline result processing tool."""


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
    help="Resfinder resistance prediction results",
)
@click.option(
    "-s",
    "--serotypefinder",
    type=click.Path(),
    help="Serotypefinder serotype prediction results",
)
@click.option("-p", "--quality", type=click.File(), help="postalignqc qc results")
@click.option("-k", "--mykrobe", type=click.File(), help="mykrobe results")
@click.option("-t", "--tbprofiler", type=click.File(), help="tbprofiler results")
@click.option(
    "--reference-genome-fasta", type=click.Path(), help="reference genome fasta file"
)
@click.option(
    "--reference-genome-gff", type=click.Path(), help="reference-genome in gff format"
)
@click.option(
    "--genome-annotation",
    type=click.Path(),
    multiple=True,
    help="Genome annotaitons bed format",
)
@click.option("--bam", type=click.Path(), help="Read mapping to reference genome")
@click.option("--snv-vcf", type=click.Path(), help="VCF with SNV variants")
@click.option("--sv-vcf", type=click.Path(), help="VCF with SV variants")
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
    serotypefinder,
    quality,
    mykrobe,
    tbprofiler,
    bam,
    reference_genome_fasta,
    reference_genome_gff,
    genome_annotation,
    snv_vcf,
    sv_vcf,
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

    if serotypefinder:
        LOG.info("Parse serotypefinder results")
        # OH typing
        res: MethodIndex | None = parse_serotypefinder_oh_typing(serotypefinder)
        if res is not None:
            results["typing_result"].extend(res)

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
                version=get_db_version(pred_res["db_version"]),
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

    # parse SNV and SV variants.
    if snv_vcf:
        results["snv_variants"] = load_variants(snv_vcf)

    if sv_vcf:
        results["sv_variants"] = load_variants(sv_vcf)

    # entries for reference genome and read mapping
    if all([bam, reference_genome_fasta, reference_genome_gff]):
        # verify that everything pertains to the same reference genome
        bam_ref_genome = get_reference_seq_accnr(bam)
        ref_accession, ref_name = get_gb_genome_version(reference_genome_gff)
        if ref_accession != bam_ref_genome:
            raise click.UsageError(
                f"Read mapping used as different reference genome; bam accnr: {bam_ref_genome}; gbff accnr: {ref_accession}"
            )

        # store file names
        fasta_idx_path = Path(f"{reference_genome_fasta}.fai")
        results["reference_genome"] = ReferenceGenome(
            name=ref_name,
            accession=ref_accession,
            fasta=Path(reference_genome_fasta).name,
            fasta_index=fasta_idx_path.name if fasta_idx_path.is_file() else None,
            genes=reference_genome_gff,
        )
        results["read_mapping"] = bam
        # add annotations
        annotations = [
            {"name": f"annotation_{i}", "file": Path(annot).name}
            for i, annot in enumerate(genome_annotation, start=1)
        ]
        for vcf in [sv_vcf, snv_vcf]:
            if vcf:
                name = "SNV" if vcf == sv_vcf else "SV"
                annotations.append({"name": name, "file": Path(vcf).name})
        # store annotation results
        results["genome_annotation"] = annotations if len(annotations) > 0 else None

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


@cli.command()
@click.option("-i", "--sample-id", required=True, help="Sample identifier")
@click.option("-b", "--bam", required=True, type=click.File(), help="bam file")
@click.option("-e", "--bed", type=click.File(), help="bed file")
@click.option("-a", "--baits", type=click.File(), help="baits file")
@click.option(
    "-r", "--reference", required=True, type=click.File(), help="reference fasta"
)
@click.option("-c", "--cpus", type=click.INT, default=1, help="cpus")
@click.option(
    "-o", "--output", required=True, type=click.File("w"), help="output filepath"
)
def create_qc_result(sample_id, bam, bed, baits, reference, cpus, output) -> None:
    """Generate QC metrics regarding bam file"""
    if bam and reference:
        LOG.info("Parse alignment results")
        parse_alignment_results(sample_id, bam, reference, cpus, output, bed, baits)
    click.secho("Finished generating QC output", fg="green")


@cli.command()
@click.option("-v", "--vcf", type=click.Path(exists=True), help="VCF file")
@click.option("-b", "--bed", type=click.Path(exists=True), help="BED file")
@click.argument("output", type=click.Path(writable=True))
def annotate_delly(vcf, bed, output):
    """Annotate Delly SV varinats with genes in BED file."""
    output = Path(output)
    # load annotation
    if bed is not None:
        annotation = pysam.TabixFile(bed, parser=pysam.asTuple())
    else:
        raise click.UsageError("You must provide a annotation file.")

    vcf_obj = VCF(vcf)
    variant = next(vcf_obj)
    annot_chrom = False
    if not variant.CHROM in annotation.contigs:
        if len(annotation.contigs) > 1:
            raise click.UsageError(
                f'"{variant.CHROM}" not in BED file and the file contains {len(annotation.contigs)} chromosomes'
            )
        # if there is only one "chromosome" in the bed file
        annot_chrom = True
        LOG.warning("Annotating variant chromosome to %s", annotation.contigs[0])
    # reset vcf file
    vcf_obj = VCF(vcf)
    vcf_obj.add_info_to_header(
        {
            "ID": "gene",
            "Description": "overlapping gene",
            "Type": "Character",
            "Number": "1",
        }
    )
    vcf_obj.add_info_to_header(
        {
            "ID": "locus_tag",
            "Description": "overlapping tbdb locus tag",
            "Type": "Character",
            "Number": "1",
        }
    )

    # open vcf writer
    writer = Writer(output.absolute(), vcf_obj)
    annotate_delly_variants(writer, vcf_obj, annotation, annot_chrom=annot_chrom)

    click.secho(f"Wrote annotated delly variants to {output}", fg="green")
