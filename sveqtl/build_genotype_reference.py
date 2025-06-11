#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""Merge SV callsets to produce a genotype reference panel."""


import os
import shutil
import tempfile
from typing import Dict, List, Tuple

import pysam  # type: ignore

from sveqtl.genotyping.catalogue_config import SOURCE_CONFIGS
from sveqtl.genotyping.ref_downloader import RefDownloader
from sveqtl.genotyping.sv_catalogue import SVCatalogue
from sveqtl.genotyping.sv_reference_merger import SVReferenceMerger


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


def _filter_vcf_ref_bounds(vcf_path: str, fai_path: str) -> Tuple[int, int]:
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
        return kept, removed


def _check_for_index(reference: str) -> str:
    """Check if the reference FASTA file has an index."""
    reference_index = f"{reference}.fai"
    if not os.path.exists(reference_index):
        pysam.faidx(reference)
        if not os.path.exists(reference_index):
            raise FileNotFoundError(
                f"Reference index {reference_index} could not be created. "
                "Ensure that the reference FASTA file is valid."
            )
    return reference_index


def build_genotype_reference(
    reference: str,
    output_path: str = ".",
    concordance: bool = False,
) -> None:
    """Produce an SV genotyping reference panel."""
    reference_index = _check_for_index(reference)

    short_read_svs: List[Dict] = []
    long_read_svs: List[Dict] = []

    # Download and preprocess SV reference callsets
    downloader = RefDownloader(configs=SOURCE_CONFIGS, output_path=output_path)
    downloader.download_ref_callsets()
    downloader.filter_callsets()

    # Normalize each VCF
    for source_name, config in SOURCE_CONFIGS.items():
        filtered_vcf_path = f"{output_path}/{config.output_name}.filtered.vcf"

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

    # Instantiate class and run the merge pipeline
    merger = SVReferenceMerger(
        short_read_svs=short_read_svs,
        long_read_svs=long_read_svs,
        ref_fasta=reference,
        concordance=concordance,
    )
    merger.merge_callsets()

    if concordance:
        merged_vcf_path = f"{output_path}/sv_genotype_reference_concordance.vcf"
    else:
        merged_vcf_path = f"{output_path}/sv_genotype_reference.vcf"

    merger.write_merged(merged_vcf_path)
    print("[RefMerger] Finished merging SV callsets.")

    # Ensure that SVs are within the reference genome bounds
    kept, removed = _filter_vcf_ref_bounds(
        vcf_path=merged_vcf_path,
        fai_path=reference_index,
    )
    print(
        "Filtered SV callset to reference genome bounds. "
        f"removed {removed} OOR variants; final count {kept}."
    )
