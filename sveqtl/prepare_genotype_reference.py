#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""Merge SV callsets to produce a genotype reference panel."""


import argparse
import os
import shutil
import tempfile
from typing import Dict, List, Tuple

import pysam

from sveqtl.catalogue_config import SOURCE_CONFIGS
from sveqtl.ref_downloader import RefDownloader
from sveqtl.sv_catalogue import SVCatalogue
from sveqtl.sv_reference_merger import SVReferenceMerger


def _load_contig_lengths(fai_path: str) -> Dict[str, int]:
    """Get contig lengths from a FASTA index.

    Returns:
        {str: int}: Dictionary mapping contig names to their lengths.
    """
    contig_lengths: Dict[str, int] = {}
    with open(fai_path, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip().split("\t")
            if len(parts) >= 2:
                contig, length = parts[0], int(parts[1])
                contig_lengths[contig] = length

    if not contig_lengths:
        raise RuntimeError(f"No entries read from {fai_path}")

    return contig_lengths


def filter_vcf_ref_bounds(vcf_path: str, fai_path: str) -> Tuple[int, int]:
    """Remove variants where pos or end are outside the reference genome.

    Returns:
      kept, removed: Number of variants kept and removed.
    """
    contig_lengths = _load_contig_lengths(fai_path)
    kept, removed = 0, 0

    with pysam.VariantFile(vcf_path) as fin:
        header = fin.header.copy()
        parent_dir = os.path.dirname(vcf_path)

        with tempfile.NamedTemporaryFile("wb", delete=False, dir=parent_dir) as tmp:
            tmp_path = tmp.name
            with pysam.VariantFile(tmp, "w", header=header) as out:
                for rec in fin:
                    chrom_len = contig_lengths.get(rec.chrom)
                    if chrom_len is not None and rec.stop <= chrom_len:
                        out.write(rec)
                        kept += 1
                    else:
                        removed += 1

        shutil.move(tmp_path, vcf_path)
        os.remove(tmp_path)
        return kept, removed


def main() -> None:
    """Produce an SV genotyping reference panel."""
    argparser = argparse.ArgumentParser(
        description="Merge SV callsets to produce a genotype reference panel."
    )
    argparser.add_argument(
        "--output_path",
        type=str,
        default=".",
        help="Path to output directory for merged SV reference panel.",
    )
    argparser.add_argument(
        "--reference_index",
        type=str,
        required=True,
        help="Path to the reference genome index file (FASTA index).",
    )
    args = argparser.parse_args()

    short_read_svs: List[Dict] = []
    long_read_svs: List[Dict] = []

    # Download and preprocess SV reference callsets
    downloader = RefDownloader(configs=SOURCE_CONFIGS, output_path=args.output_path)
    downloader.download_ref_callsets()
    downloader.filter_callsets()

    # Normalize each VCF
    for source_name, config in SOURCE_CONFIGS.items():
        filtered_vcf_path = f"{args.output_path}/{config.output_name}.filtered.vcf"

        catalogue = SVCatalogue(
            vcf_path=filtered_vcf_path,
            source_name=source_name,
            priority=config.priority,
            read_type=config.read_type,
            source_rules={source_name: config.rules.__dict__},
        )
        if catalogue.read_type == "short_read":
            short_read_svs.extend(catalogue.variants)
        else:
            long_read_svs.extend(catalogue.variants)

    # Instantiate the SVReferenceMerger class and run the merge pipeline
    merger = SVReferenceMerger(
        short_read_svs=short_read_svs, long_read_svs=long_read_svs
    )
    merger.merge_callsets()

    merged_vcf_path = f"{args.output_path}/sv_genotype_reference.vcf"
    merger.write_merged(merged_vcf_path)
    print("[RefMerger] Finished merging SV callsets.")

    # Filter step to ensure that all SVs are within the reference genome bounds
    kept, removed = filter_vcf_ref_bounds(
        vcf_path=merged_vcf_path,
        fai_path=args.reference_index,
    )
    print(
        "Filtered SV callset to reference genome bounds. "
        f"removed {removed} OOR variants; final count {kept}."
    )


if __name__ == "__main__":
    main()
