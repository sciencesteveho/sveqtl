"""Command-line interface for sveqtl analysis."""

import argparse
import logging
from pathlib import Path
import subprocess
import sys

from sveqtl.build_genotype_reference import build_genotype_reference
from sveqtl.filter_genotypes import combine_and_filter_genotypes
from sveqtl.utils.parse_sample_manifests import parse_sample_manifests
from sveqtl.utils.parsers import _build_genotype_parser
from sveqtl.utils.parsers import _build_genotype_reference_parser
from sveqtl.utils.parsers import _build_idxdepth_parser
from sveqtl.utils.parsers import _filter_genotype_parser

_LOG_FORMAT = "%(asctime)s | %(levelname)-8s | " "%(name)s:%(lineno)d - %(message)s"
logging.basicConfig(
    level=logging.INFO,
    format=_LOG_FORMAT,
    datefmt="%H:%M:%S",
)
LOGGER = logging.getLogger("sveqtl.cli")


def _build_genotype_reference(args: argparse.Namespace) -> None:
    """Entry point for building the SV genotyping reference."""
    LOGGER.info("Building SV genotyping reference panel.")
    build_genotype_reference(
        reference=args.reference,
        output_path=args.output_dir,
    )
    LOGGER.info(f"SV genotyping reference panel built and saved to {args.output_dir}.")


def _run_idx_depth(args: argparse.Namespace) -> None:
    """Entry point for running idxdepth to prepare sample manifests."""
    script = Path(__file__).parent.parent / "shell" / "run_idxdepth.sh"
    cmd = [
        "bash",
        str(script),
        args.paragraph_dir,
        args.bam_file,
        args.reference,
        str(args.threads),
        args.idxdepth_out,
    ]
    LOGGER.info(f"Running idx depth via: {' '.join(cmd)}")

    result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, text=True)
    sample_id = result.stdout.strip()
    LOGGER.info(f"Idxdepth complete for sample: {sample_id}")

    LOGGER.info("Parsing sample manifest JSON → TSV")
    parse_sample_manifests(
        idxdepth_out=args.idxdepth_out,
        sample=sample_id,
    )
    LOGGER.info(
        f"Sample manifest written to: {args.idxdepth_out}/manifests/{sample_id}.manifest.tsv"
    )


def _genotype(args: argparse.Namespace) -> None:
    """Entry point for genotyping SVs with Paragraph."""
    script = Path(__file__).parent.parent / "shell" / "paragraph_genotyping.sh"
    LOGGER.info(f"Invoking Paragraph genotyping script: {script}")
    subprocess.check_call(
        [
            "bash",
            str(script),
            args.sample,
            args.paragraph_dir,
            args.manifest_dir,
            args.reference,
            args.genotype_out_dir,
            str(args.threads),
        ]
    )
    LOGGER.info("Paragraph genotyping complete.")


def _combine_and_filter_genotypes(args: argparse.Namespace) -> None:
    """Entry point for filtering genotyped SVs."""
    LOGGER.info("Combining and filtering SV VCFs for QTL analysis.")
    combine_and_filter_genotypes(
        sample_dir=Path(args.sample_dir),
        threads=args.threads,
        min_call_rate=args.min_call_rate,
        hwe_alpha=args.hwe_alpha,
        min_maf=args.min_maf,
    )


def main() -> None:
    """Main entry point for sveqtl analysis CLI."""
    parser = argparse.ArgumentParser(
        prog="sveqtl",
        description="SV-eQTL analysis toolkit",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    _build_genotype_reference_parser(subparsers, _build_genotype_reference)
    _build_idxdepth_parser(subparsers, _run_idx_depth)
    _build_genotype_parser(subparsers, _genotype)
    _filter_genotype_parser(subparsers, _combine_and_filter_genotypes)

    args = parser.parse_args()

    try:
        args.func(args)
    except Exception as exc:
        LOGGER.exception(f"Unhandled exception: {exc}")
        LOGGER.error("Exiting with error.")
        sys.exit(1)


if __name__ == "__main__":
    main()
