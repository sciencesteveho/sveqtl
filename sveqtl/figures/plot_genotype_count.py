#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Plot # of genotyped SVs per sample in cohort."""

import argparse
from collections import Counter
from collections import defaultdict
from typing import Dict, List, Optional

from matplotlib import cm
import matplotlib.pyplot as plt
import pandas as pd  # type: ignore
import pysam

from sveqtl.figures import _set_matplotlib_publication_parameters


def count_sv_genotypes(
    vcf_path: str,
    sv_types: List[str],
) -> Dict[str, Dict[str, int]]:
    """Scan a multi-sample VCF and count, for each sample, how many times
    it carries an SV of each type.
    """
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)

    # initialize counts per sample
    counts = {s: {sv: 0 for sv in sv_types} for s in samples}

    for rec in vcf:
        sv = rec.info.get("SVTYPE")
        if sv not in sv_types:
            continue
        for s in samples:
            gt = rec.samples[s]["GT"]
            if gt is not None and any(a == 1 for a in gt if a is not None):
                counts[s][sv] += 1
    return counts


def count_unique_sv_loci(
    vcf_path: str,
    sv_types: List[str],
) -> int:
    """Count the number of unique SV loci with at least one non-reference
    call.
    """
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)
    unique_count = 0

    for rec in vcf:
        sv = rec.info.get("SVTYPE")
        if sv not in sv_types:
            continue

        for s in samples:
            gt = rec.samples[s]["GT"]
            if gt is not None and any(a == 1 for a in gt if a is not None):
                unique_count += 1
                break
    return unique_count


def build_sv_dataframe(counts: Dict[str, Dict[str, int]]) -> pd.DataFrame:
    """Turn the nested count dict into a DataFrame, add total_svs, and sort."""
    df = pd.DataFrame.from_dict(counts, orient="index")
    df["total_svs"] = df.sum(axis=1)
    df = df.sort_values("total_svs").reset_index().rename(columns={"index": "sample"})
    return df


def plot_sv_counts(
    df: pd.DataFrame,
    sv_types: List[str],
    figsize: tuple = (10, 6),
    cmap_name: str = "viridis",
) -> None:
    """Scatter-plot per-sample SV counts by type ordered by total_svs."""
    if sv_types is None:
        sv_types = [c for c in df.columns if c not in ("sample", "total_svs")]

    x = df.index
    cmap = cm.get_cmap(cmap_name, len(sv_types))

    plt.figure(figsize=figsize)
    for i, sv in enumerate(sv_types):
        plt.scatter(
            x,
            df[sv],
            label=sv,
            color=cmap(i),
            s=50,
            alpha=0.8,
        )
    plt.xlabel("Sample (sorted by total SV count)")
    plt.ylabel("Count of structural variants")
    plt.title("Per-sample SV counts by type")
    plt.legend(title="SV type")
    plt.tight_layout()
    plt.savefig("sv_counts_per_sample.svg")


def main(vcf: str = "cohort.genotypes_filtered.vcf") -> None:
    """Main function to parse arguments, print summary, and plot results."""

    _set_matplotlib_publication_parameters()
    sv_types = ["DEL", "DUP", "INS", "INV"]

    counts = count_sv_genotypes(
        vcf_path=vcf,
        sv_types=sv_types,
    )
    df = build_sv_dataframe(counts)

    total_sv_calls = int(df["total_svs"].sum())
    unique_sv_loci = count_unique_sv_loci(
        vcf_path=vcf,
        sv_types=sv_types,
    )
    print(f"Total SV calls across samples: {total_sv_calls}")
    print(f"Total unique SV loci with at least one call: {unique_sv_loci}")

    plot_sv_counts(
        df=df,
        sv_types=sv_types,
    )


if __name__ == "__main__":
    main()
