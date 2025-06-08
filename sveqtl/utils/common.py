#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Utilities for sv-eqtl methods and analysis."""

from pathlib import Path
import subprocess
from typing import List


def _glob_files(sample_dir: Path, pattern: str) -> List[Path]:
    """Find all files matching pattern in the given directory."""
    return [file for file in sample_dir.glob(pattern) if file.is_file()]


def merge_vcfs(vcfs: List[Path], out_vcf_gz: Path, threads: int) -> None:
    """Merge multiple VCFs into a single cohort VCF."""
    if out_vcf_gz.exists():
        return

    subprocess.check_call(
        [
            "bcftools",
            "merge",
            "--threads",
            str(threads),
            "-Oz",
            "-o",
            str(out_vcf_gz),
            *map(str, vcfs),
        ]
    )
    subprocess.check_call(["bcftools", "index", "-t", str(out_vcf_gz)])
