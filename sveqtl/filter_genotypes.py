#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Filter and merge SV genotypes."""


from dataclasses import dataclass
from pathlib import Path
import subprocess
from typing import Dict, List

from cyvcf2 import VCF  # type: ignore
import pysam  # type: ignore
from snphwe import snphwe  # type: ignore

from sveqtl.utils.common import _glob_files
from sveqtl.utils.common import _index_vcf


@dataclass(frozen=True)
class GenotypeFilters:
    """Configuration class for filtering thresholds."""

    min_call_rate: float = 0.50
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
    def call_rate(self) -> float:
        """Fraction of called (non-missing) genotypes."""
        total = self.total_called + self.n_missing
        return self.total_called / total if total else 0.0

    @property
    def maf(self) -> float:
        """Minor allele frequency."""
        if self.total_called == 0:
            return 0.0
        alt_count = self.n_het + 2 * self.n_alt
        alt_af = alt_count / (2 * self.total_called)
        return min(alt_af, 1.0 - alt_af)


def combine_and_filter_genotypes(
    sample_dir: Path,
    threads: int = 16,
    min_call_rate: float = 0.50,
    hwe_alpha: float = 1e-4,
    min_maf: float = 0.05,
    merged_vcf: str = "cohort.merged.vcf.gz",
    out_with_maf: str = "cohort.genotypes_filtered.maf.vcf",
    out_without_maf: str = "cohort.genotypes_filtered.vcf",
) -> None:
    """Entry point for genotype QC pipeline:
      1. Locate all per-sample VCFs under --sample-dir.
      2. Merge the cleaned VCFs into a single cohort VCF.
      3. Filter the cohort VCF for call rate, MAF, and HWE.

    Examples::
      >>> combine_and_filter_genotypes(
      ...   sample_dir=Path(args.sample_dir),
      ...   threads=args.threads,
      ...   min_call_rate=args.min_call_rate,
      ...   hwe_alpha=args.hwe_alpha,
      ...   min_maf=args.min_maf,
      >>>  )
    """
    genotype_filters = GenotypeFilters(
        min_call_rate=min_call_rate,
        hwe_alpha=hwe_alpha,
        min_maf=min_maf,
    )

    vcfs = _glob_files(
        sample_dir=sample_dir,
        pattern="*/genotypes.vcf.gz",
    )
    if not vcfs:
        raise SystemExit("No per-sample VCFs found. Check --sample_dir path")

    for vcf in vcfs:
        _index_vcf(vcf, threads=threads)

    print(f"Found {len(vcfs)} per-sample VCFs")

    _merge_vcfs(
        vcfs=vcfs,
        out_vcf_gz=merged_vcf,
        threads=threads,
    )
    print(f"Merged cohort VCF written to: {merged_vcf}")

    _filter_population_vcf(
        population_vcf=merged_vcf,
        out_with_maf=out_with_maf,
        out_without_maf=out_without_maf,
        genotype_filters=genotype_filters,
        threads=threads,
    )
    print(f"PASS VCF with MAF filter written to: {out_with_maf}")
    print(f"PASS VCF without MAF filter written to: {out_without_maf}")


def _test_hardy_weinberg_equilibrium(counts: GenotypeCounts) -> float:
    """Calculate HWE p-value from genotype counts."""
    if counts.total_called == 0:
        return 1.0

    return snphwe(counts.n_het, counts.n_alt, counts.n_ref)


def _merge_vcfs(
    vcfs: List[Path],
    out_vcf_gz: str,
    threads: int,
) -> None:
    """Merge multiple VCFs into a single cohort VCF."""
    merge_cmd = [
        "bcftools",
        "merge",
        "--threads",
        str(threads),
        "-Oz",
        *map(str, vcfs),
    ]

    merge_proc = subprocess.Popen(merge_cmd, stdout=subprocess.PIPE)
    with open(out_vcf_gz, "wb") as f:
        f.write(merge_proc.communicate()[0])

    if merge_ret := merge_proc.wait():
        raise subprocess.CalledProcessError(merge_ret, merge_cmd)

    subprocess.check_call(["bcftools", "index", "-t", out_vcf_gz])


def _count_variant_genotypes(variant) -> GenotypeCounts:
    """Count genotypes for a variant using cyvcf2 built-in properties."""
    return GenotypeCounts(
        n_ref=variant.num_hom_ref,
        n_het=variant.num_het,
        n_alt=variant.num_hom_alt,
        n_missing=variant.num_unknown,
    )


def _keep_only_gt(rec: pysam.VariantRecord) -> None:
    """Strip out all sample fields except GT in-place."""
    for _, record in rec.samples.items():
        gt = record["GT"]
        record.clear()
        record["GT"] = gt


def _filter_population_vcf(
    population_vcf: str,
    out_with_maf: str,
    out_without_maf: str,
    genotype_filters: GenotypeFilters,
    threads: int,
) -> None:
    """Filter merged, sample-level VCFs by call rate, HWE, and MAF."""

    vcf_reader = VCF(population_vcf)
    with pysam.VariantFile(population_vcf, threads=threads) as ivcf:
        header = ivcf.header.copy()

        with pysam.VariantFile(
            out_with_maf, "w", header=header, threads=threads  # type: ignore
        ) as o_maf, pysam.VariantFile(
            out_without_maf, "w", header=header, threads=threads  # type: ignore
        ) as o_no_maf:

            stats: Dict[str, Dict[str, int]] = {
                "with": {"kept": 0, "dropped": 0},
                "without": {"kept": 0, "dropped": 0},
            }

            dropped_call_rate = dropped_hwe = dropped_maf = 0

            pysam_records = ivcf.fetch()
            for cyvcf_variant in vcf_reader:
                try:
                    pysam_rec = next(pysam_records)
                except StopIteration:
                    break

                counts = _count_variant_genotypes(cyvcf_variant)
                call_rate_fail = counts.call_rate < genotype_filters.min_call_rate
                hwe_fail = (
                    _test_hardy_weinberg_equilibrium(counts)
                    < genotype_filters.hwe_alpha
                )
                maf_fail = counts.maf < genotype_filters.min_maf

                if call_rate_fail:
                    dropped_call_rate += 1
                if hwe_fail:
                    dropped_hwe += 1
                if maf_fail:
                    dropped_maf += 1

                keep_maf = not (call_rate_fail or hwe_fail or maf_fail)
                keep_no_maf = not (call_rate_fail or hwe_fail)

                _keep_only_gt(pysam_rec)
                if keep_maf:
                    o_maf.write(pysam_rec)
                    stats["with"]["kept"] += 1
                else:
                    stats["with"]["dropped"] += 1

                if keep_no_maf:
                    o_no_maf.write(pysam_rec)
                    stats["without"]["kept"] += 1
                else:
                    stats["without"]["dropped"] += 1

    _print_filtering_stats(
        dropped_call_rate=dropped_call_rate,
        dropped_hwe=dropped_hwe,
        dropped_maf=dropped_maf,
        stats=stats,
    )


def _print_filtering_stats(
    dropped_call_rate: int,
    dropped_hwe: int,
    dropped_maf: int,
    stats: Dict[str, Dict[str, int]],
) -> None:
    """Print filtering statistics."""
    print(
        f"{dropped_call_rate:,} removed due to call rate, "
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
