"""Argument parsers for sveqtl CLI commands."""

import argparse
from typing import Callable


def _build_genotype_reference_parser(
    subparsers: argparse._SubParsersAction,
    func: Callable,
) -> argparse.ArgumentParser:
    """Parser for the build_genotype_ref command."""
    build_ref_parser = subparsers.add_parser(
        "build_genotype_ref",
        help="Build the SV genotyping reference panel.",
    )
    build_ref_parser.add_argument(
        "--reference",
        required=True,
        help="Path to reference FASTA. Ensure that an index is in the same directory.",
    )
    build_ref_parser.add_argument(
        "--output_dir",
        required=True,
        help="Output location for the SV genotyping reference panel.",
    )
    build_ref_parser.add_argument(
        "--concordance",
        action="store_true",
        help="Only keep variants supported by 2+ datasets.",
    )
    build_ref_parser.set_defaults(func=func)
    return build_ref_parser


def _build_idxdepth_parser(
    subparsers: argparse._SubParsersAction,
    func: Callable,
) -> argparse.ArgumentParser:
    """Parser for the run_idx_depth command."""
    idxdepth_parser = subparsers.add_parser(
        "run_idx_depth",
        help="Run index-depth analysis (wraps shell/run_idxdepth.sh)",
    )
    idxdepth_parser.add_argument(
        "--paragraph_dir",
        type=str,
        required=True,
        help="Path to Paragraph installation.",
    )
    idxdepth_parser.add_argument(
        "--bam_file", required=True, help="CRAM/BAM file to process."
    )
    idxdepth_parser.add_argument(
        "--reference",
        required=True,
        help="Reference FASTA.",
    )
    idxdepth_parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads to use.",
    )
    idxdepth_parser.add_argument(
        "--idxdepth_out",
        required=True,
        help="Output directory for idxdepth results.",
    )
    idxdepth_parser.set_defaults(func=func)
    return idxdepth_parser


def _build_genotype_parser(
    subparsers: argparse._SubParsersAction,
    func: Callable,
) -> argparse.ArgumentParser:
    """Parser for the genotype command."""
    geno_parser = subparsers.add_parser(
        "genotype",
        help="Run Paragraph genotyping (wraps shell/paragraph_genotyping.sh)",
    )
    geno_parser.add_argument(
        "--paragraph_dir",
        required=True,
        help="Path to the Paragraph installation directory.",
    )
    geno_parser.add_argument(
        "--manifest_dir",
        required=True,
        help="Directory containing per-sample manifest TSVs.",
    )
    geno_parser.add_argument("--reference", required=True, help="Reference FASTA.")
    geno_parser.add_argument(
        "--threads", type=int, default=32, help="Number of threads to use."
    )
    geno_parser.add_argument(
        "--genotype_out_dir",
        required=True,
        help="Output directory for pilot genotyping.",
    )
    geno_parser.add_argument(
        "sample", help="Sample name (must match manifest basename)."
    )
    geno_parser.set_defaults(func=func)
    return geno_parser


def _filter_genotype_parser(
    subparsers: argparse._SubParsersAction,
    func: Callable,
) -> argparse.ArgumentParser:
    """Parser for the combine_and_filter_genotypes command."""
    combine_parser = subparsers.add_parser(
        "combine_and_filter_genotypes",
        help="Merge & filter genotype VCFs.",
    )
    combine_parser.add_argument(
        "--sample_dir",
        required=True,
        help="Directory of per-sample VCFs.",
    )
    combine_parser.add_argument(
        "--threads",
        type=int,
        default=18,
        help="Number of threads to use.",
    )
    combine_parser.add_argument(
        "--min_call_rate",
        type=float,
        default=0.50,
        help="Maximum allowed missing genotype fraction.",
    )
    combine_parser.add_argument(
        "--hwe_alpha",
        type=float,
        default=1e-4,
        help="HWE significance threshold.",
    )
    combine_parser.add_argument(
        "--min_maf",
        type=float,
        default=0.05,
        help="Minimum minor allele frequency.",
    )
    combine_parser.set_defaults(func=func)
    return combine_parser
