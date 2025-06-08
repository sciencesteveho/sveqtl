#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Filter and merge SV genotypes."""


import argparse
from dataclasses import dataclass
from pathlib import Path
import subprocess
from typing import Dict, Iterable, List, Tuple

import pysam  # type: ignore
from snphwe import snphwe  # type: ignore

from sveqtl.utils.common import _glob_files
from sveqtl.utils.common import _index_vcf


@dataclass(frozen=True)
class GenotypeFilters:
    """Configuration class for filtering thresholds."""

    max_missing_frac: float = 0.50
    hwe_alpha: float = 1e-4
    min_maf: float = 0.05


@dataclass
class GenotypeCounts:
    """Counts of genotypes for a single variant."""

    n_ref: int
    n_het: int
    n_alt: int
    n_missing: int

    @property
    def total_called(self) -> int:
        """Total number of called (not-missing) genotypes."""
        return self.n_ref + self.n_het + self.n_alt

    @property
    def missing_fraction(self) -> float:
        """Fraction of missing genotypes."""
        total = self.total_called + self.n_missing
        return self.n_missing / total if total else 1.0

    @property
    def maf(self) -> float:
        """Minor allele frequency."""
        if self.total_called == 0:
            return 0.0
        alt_count = self.n_het + 2 * self.n_alt
        alt_af = alt_count / (2 * self.total_called)
        return min(alt_af, 1.0 - alt_af)


def _test_hardy_weinberg_equilibrium(counts: GenotypeCounts) -> float:
    """Calculate HWE p-value from genotype counts."""
    return snphwe(counts.n_het, counts.n_alt, counts.n_ref)


def merge_vcfs(
    vcfs: List[Path],
    out_vcf_gz: Path,
    threads: int,
) -> None:
    """Merge multiple VCFs into a single cohort VCF."""
    if out_vcf_gz.exists():
        return

    merge_cmd = [
        "bcftools",
        "merge",
        "--threads",
        str(threads),
        "-Ov",
        *map(str, vcfs),
    ]
    compress_cmd = ["bgzip", "-c", "-@", str(threads)]

    # bcftools merge | bgzip
    merge_proc = subprocess.Popen(merge_cmd, stdout=subprocess.PIPE)
    compress_proc = subprocess.Popen(
        compress_cmd, stdin=merge_proc.stdout, stdout=subprocess.PIPE
    )
    merge_proc.stdout.close()  # type: ignore

    with open(out_vcf_gz, "wb") as f:
        f.write(compress_proc.communicate()[0])

    if merge_ret := merge_proc.wait():
        raise subprocess.CalledProcessError(merge_ret, merge_cmd)
    if compress_ret := compress_proc.wait():
        raise subprocess.CalledProcessError(compress_ret, compress_cmd)

    subprocess.check_call(["bcftools", "index", "-t", str(out_vcf_gz)])


def _count_variant_genotypes(rec: pysam.VariantRecord) -> GenotypeCounts:
    """Count genotypes for a variant record."""
    n_ref = n_het = n_alt = n_miss = 0
    for sm in rec.samples.values():
        gt = sm["GT"]

        # missing genotypes
        if None in gt or -1 in gt:
            n_miss += 1
            continue

        if len(gt) == 1:  # haploid
            allele = gt[0]
            if allele == 0:
                n_ref += 1
            elif allele == 1:
                n_alt += 1
            else:
                n_miss += 1

        elif len(gt) == 2:  # diploid
            g1, g2 = gt
            if g1 == g2 == 0:
                n_ref += 1
            elif g1 == g2 == 1:
                n_alt += 1
            else:
                n_het += 1
        else:
            n_miss += 1

    return GenotypeCounts(n_ref, n_het, n_alt, n_miss)


def _prune_to_gt(header: pysam.VariantHeader) -> pysam.VariantHeader:
    """Return a copy of vcf header with all FORMAT fields except GT removed."""
    new_header = pysam.VariantHeader()
    for rec in header.records:
        if rec.type != "FORMAT" or rec["ID"] == "GT":
            new_header.add_record(rec)

    for sample in header.samples:
        new_header.add_sample(sample)

    return new_header


def _keep_only_gt(rec: pysam.VariantRecord) -> None:
    """Strip out all sample fields except GT in-place."""
    for _, record in rec.samples.items():
        gt = record["GT"]
        record.clear()
        record["GT"] = gt


def filter_population_vcf(
    population_vcf: Path,
    out_with_maf: Path,
    out_without_maf: Path,
    genotype_filters: GenotypeFilters,
    threads: int,
) -> None:
    """Filter merged, sample-level VCFs by missingness, HWE, and MAF."""
    with pysam.VariantFile(str(population_vcf), threads=threads) as ivcf:
        pruned_header = _prune_to_gt(ivcf.header)

        with pysam.VariantFile(
            str(out_with_maf), "w", header=pruned_header, threads=threads
        ) as o_maf, pysam.VariantFile(
            str(out_without_maf), "w", header=pruned_header, threads=threads
        ) as o_no_maf:

            stats: Dict[str, Dict[str, int]] = {
                "with": {"kept": 0, "dropped": 0},
                "without": {"kept": 0, "dropped": 0},
            }

            dropped_missing = dropped_hwe = dropped_maf = 0
            for rec in ivcf.fetch():
                counts = _count_variant_genotypes(rec)

                missing_fail = (
                    counts.missing_fraction > genotype_filters.max_missing_frac
                )
                hwe_fail = (
                    _test_hardy_weinberg_equilibrium(counts)
                    < genotype_filters.hwe_alpha
                )
                maf_fail = counts.maf < genotype_filters.min_maf

                if missing_fail:
                    dropped_missing += 1
                if hwe_fail:
                    dropped_hwe += 1
                if maf_fail:
                    dropped_maf += 1

                keep_maf = not (missing_fail or hwe_fail or maf_fail)
                keep_no_maf = not (missing_fail or hwe_fail)

                _keep_only_gt(rec)

                if keep_maf:
                    o_maf.write(rec)
                    stats["with"]["kept"] += 1
                else:
                    stats["with"]["dropped"] += 1

                if keep_no_maf:
                    o_no_maf.write(rec)
                    stats["without"]["kept"] += 1
                else:
                    stats["without"]["dropped"] += 1

    _print_filtering_stats(
        dropped_missing=dropped_missing,
        dropped_hwe=dropped_hwe,
        dropped_maf=dropped_maf,
        stats=stats,
    )


def _print_filtering_stats(
    dropped_missing: int,
    dropped_hwe: int,
    dropped_maf: int,
    stats: Dict[str, Dict[str, int]],
) -> None:
    """Print filtering statistics."""
    print(
        f"{dropped_missing:,} removed due to missing fraction, "
        f"{dropped_hwe:,} dropped for violating HWE, "
        f"{dropped_maf:,} for MAF threshold."
    )
    for key in ("with", "without"):
        kept = stats[key]["kept"]
        dropped = stats[key]["dropped"]
        pct = (dropped / (kept + dropped) * 100) if (kept + dropped) else 0.0
        print(
            f"{'With' if key == 'with' else 'Without'} MAF filter → "
            f"kept {kept:,}, dropped {dropped:,} ({pct:.1f}%)"
        )


def build_parser() -> argparse.ArgumentParser:
    """Build CLI parser."""
    parser = argparse.ArgumentParser("SV genotyping QC pipeline.")
    parser.add_argument(
        "--sample_dir",
        required=True,
        type=Path,
        help="Root directory containing paragraph output (one folder per sample)",
    )
    parser.add_argument(
        "--threads", default=16, type=int, help="Number of HTSlib threads"
    )
    parser.add_argument(
        "--out_with_maf",
        type=Path,
        default=Path("cohort.genotypes_filtered.maf.vcf.gz"),
        help="Output bgzipped VCF with MAF filtering",
    )
    parser.add_argument(
        "--out_without_maf",
        type=Path,
        default=Path("cohort.genotypes_filtered.vcf.gz"),
        help="Output bgzipped VCF without MAF filtering",
    )
    parser.add_argument(
        "--merged_vcf",
        type=Path,
        default=Path("cohort.merged.vcf.gz"),
        help="Intermediate merged (unfiltered) VCF",
    )
    parser.add_argument(
        "--max_missing_frac",
        type=float,
        default=0.50,
        help="Maximum missing genotype fraction per variant",
    )
    parser.add_argument(
        "--hwe_alpha",
        type=float,
        default=1e-4,
        help="One-sided HWE p-value threshold",
    )
    parser.add_argument(
        "--min_maf",
        type=float,
        default=0.05,
        help="Minimum minor-allele frequency threshold",
    )
    return parser


def main() -> None:
    """Entry point for genotype QC pipeline and performs the following:
      1. Locate all per-sample VCFs under `--sample-dir`.
      2. Merge the cleaned VCFs into a single cohort VCF and filter out
         variants that fail QC.
      3. Filter the cohort VCF for genotyping fraction, MAF, and HWE.
      
    Examples::
      >>> python filter_genotypes.py \
      ...     --sample-dir /path/to/samples \
      ...     --threads 16
    """
    args = build_parser().parse_args()

    genotype_filters = GenotypeFilters(
        max_missing_frac=args.max_missing_frac,
        hwe_alpha=args.hwe_alpha,
        min_maf=args.min_maf,
    )

    vcfs = _glob_files(
        sample_dir=args.sample_dir,
        pattern="*/genotypes_filtered.vcf.gz",
    )
    if not vcfs:
        raise SystemExit("No per-sample VCFs found. Check --sample_dir path")

    for vcf in vcfs:
        _index_vcf(vcf, threads=args.threads)

    print(f"Found {len(vcfs)} per-sample VCFs")

    merge_vcfs(
        vcfs=vcfs,
        out_vcf_gz=args.merged_vcf,
        threads=args.threads,
    )
    print(f"Merged cohort VCF written to: {args.merged_vcf}")

    filter_population_vcf(
        population_vcf=args.merged_vcf,
        out_with_maf=args.out_with_maf,
        out_without_maf=args.out_without_maf,
        genotype_filters=genotype_filters,
        threads=args.threads,
    )
    print(f"PASS VCF with MAF filter written to: {args.out_with_maf}")
    print(f"PASS VCF without MAF filter written to: {args.out_without_maf}")


if __name__ == "__main__":
    main()
