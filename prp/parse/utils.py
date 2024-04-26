"""Shared utility functions."""
import os
import io
from datetime import datetime
from typing import Tuple

from ..models.phenotype import (
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


def reformat_date_str(input_date: str) -> str:
    """Reformat date string into DDMMYY format"""
    # Parse the date string
    try:
        parsed_date = datetime.strptime(input_date, "%Y-%m-%d %H:%M:%S.%f")
    except ValueError:
        parsed_date = datetime.strptime(input_date, "%a %b %d %H:%M:%S %Y %z")

    # Format as DDMMYY
    formatted_date = parsed_date.date().isoformat()
    return formatted_date


def get_db_version(db_version: dict) -> str:
    """Get database version"""
    backup_version = db_version["name"] + "_" + reformat_date_str(db_version["Date"])
    return db_version["commit"] if "commit" in db_version else backup_version


def _get_path(symlink_dir: str, subdir: str, filepath: str) -> str:
    """Get absolute/symlink path"""
    return os.path.join(symlink_dir, subdir, filepath) if symlink_dir else os.path.realpath(filepath)


def parse_input_dir(input_dir: str, jasen_dir: str, output_dir: str):
    input_arrays = []
    input_dir = input_dir.rstrip("/")
    species = input_dir.split("/")[-1]
    softwares = os.listdir(input_dir)
    output_dir = os.path.join(input_dir, "analysis_result") if not output_dir else output_dir
    if os.path.exists(input_dir):
        analysis_results_dir = os.path.join(input_dir, "analysis_result")
        for filename in os.listdir(analysis_results_dir):
            if filename.endswith(".json"):
                sample_id = filename.rstrip("_result.json")
                sample_array = create_sample_array(species, input_dir, jasen_dir, sample_id, output_dir)
                input_arrays.append(sample_array)
    return input_arrays


def create_sample_array(species, input_dir, jasen_dir, sample_id, output_dir):
    output = os.path.abspath(os.path.join(output_dir, f"{sample_id}_result.json"))
    bam = os.path.abspath(os.path.join(input_dir, f"bam/{sample_id}.bam"))
    kraken = os.path.abspath(os.path.join(input_dir, f"kraken/{sample_id}_bracken.out"))
    quality = os.path.abspath(os.path.join(input_dir, f"postalignqc/{sample_id}_bwa.qc"))
    quast = os.path.abspath(os.path.join(input_dir, f"quast/{sample_id}_quast.tsv"))
    run_metadata = os.path.abspath(os.path.join(input_dir, f"analysis_metadata/{sample_id}_analysis_meta.json"))
    if species == "mtuberculosis":
        reference_genome_fasta = os.path.abspath(os.path.join(jasen_dir, "assets/genomes/mycobacterium_tuberculosis/NC_000962.3.fasta"))
        reference_genome_gff = os.path.abspath(os.path.join(jasen_dir, "assets/genomes/mycobacterium_tuberculosis/NC_000962.3.gff"))
        sv_vcf = os.path.abspath(os.path.join(input_dir, f"annotate_delly/{sample_id}_annotated_delly.vcf"))
        mykrobe = os.path.abspath(os.path.join(input_dir, f"mykrobe/{sample_id}_mykrobe.csv"))
        tbprofiler = os.path.abspath(os.path.join(input_dir, f"tbprofiler_mergedb/{sample_id}_tbprofiler.json"))
        return {
            "output": output,
            "bam": bam,
            "kraken": kraken,
            "quality": quality,
            "quast": quast,
            "reference_genome_fasta": reference_genome_fasta,
            "reference_genome_gff": reference_genome_gff,
            "run_metadata": run_metadata,
            "sv_vcf": sv_vcf,
            "mykrobe": mykrobe,
            "tbprofiler": tbprofiler,
            "amrfinder": None,
            "cgmlst": None,
            "mlst": None,
            "resfinder": None,
            "serotypefinder": None,
            "virulencefinder": None,
            "process_metadata": None,
        }
    elif species == "saureus" or species == "ecoli" or species == "kpneumoniae":
        process_metadata = []
        amrfinder = os.path.abspath(os.path.join(input_dir, f"amrfinderplus/{sample_id}_amrfinder.out"))
        cgmlst = os.path.abspath(os.path.join(input_dir, f"chewbbaca/{sample_id}_chewbbaca.out"))
        mlst = os.path.abspath(os.path.join(input_dir, f"mlst/{sample_id}_mlst.json"))
        resfinder = os.path.abspath(os.path.join(input_dir, f"resfinder/{sample_id}_resfinder.json"))
        resfinder_meta = os.path.abspath(os.path.join(input_dir, f"resfinder/{sample_id}_resfinder_meta.json"))
        serotypefinder = os.path.abspath(os.path.join(input_dir, f"serotypefinder/{sample_id}_serotypefinder.json"))
        serotypefinder_meta = os.path.abspath(os.path.join(input_dir, f"serotypefinder/{sample_id}_serotypefinder_meta.json"))
        virulencefinder = os.path.abspath(os.path.join(input_dir, f"virulencefinder/{sample_id}_virulencefinder.json"))
        virulencefinder_meta = os.path.abspath(os.path.join(input_dir, f"virulencefinder/{sample_id}_virulencefinder_meta.json"))
        process_metadata.append(resfinder_meta)
        process_metadata.append(serotypefinder_meta)
        process_metadata.append(virulencefinder_meta)
        if species == "saureus":
            reference_genome_fasta = os.path.abspath(os.path.join(jasen_dir, "assets/genomes/staphylococcus_aureus/NC_002951.2.fasta"))
            reference_genome_gff = os.path.abspath(os.path.join(jasen_dir, "assets/genomes/staphylococcus_aureus/NC_002951.2.gff"))
        if species == "ecoli":
            reference_genome_fasta = os.path.abspath(os.path.join(jasen_dir, "assets/genomes/escherichia_coli/NC_000913.3.fasta"))
            reference_genome_gff = os.path.abspath(os.path.join(jasen_dir, "assets/genomes/escherichia_coli/NC_000913.3.gff"))
        if species == "kpneumoniae":
            reference_genome_fasta = os.path.abspath(os.path.join(jasen_dir, "assets/genomes/klebsiella_pneumoniae/NC_016845.1.fasta"))
            reference_genome_gff = os.path.abspath(os.path.join(jasen_dir, "assets/genomes/klebsiella_pneumoniae/NC_016845.1.gff"))
        return {
            "output": output,
            "bam": bam,
            "kraken": kraken,
            "quality": quality,
            "quast": quast,
            "reference_genome_fasta": reference_genome_fasta,
            "reference_genome_gff": reference_genome_gff,
            "run_metadata": run_metadata,
            "sv_vcf": None,
            "mykrobe": None,
            "tbprofiler": None,
            "amrfinder": amrfinder,
            "cgmlst": cgmlst,
            "mlst": mlst,
            "resfinder": resfinder,
            "serotypefinder": serotypefinder,
            "virulencefinder": virulencefinder,
            "process_metadata": process_metadata,
        }
