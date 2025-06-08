#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Utilities for sv-eqtl methods and analysis."""

from pathlib import Path
import subprocess
from typing import List


def _glob_files(sample_dir: Path, pattern: str) -> List[Path]:
    """Find all files matching pattern in the given directory."""
    return [file for file in sample_dir.glob(pattern) if file.is_file()]


def _index_vcf(vcf: Path, threads: int = 1) -> None:
    """Index a bgzipped VCF if it does not already have an index."""
    tbi = Path(f"{str(vcf)}.tbi")
    if not tbi.exists():
        subprocess.check_call(
            [
                "bcftools",
                "index",
                "-f",
                "--threads",
                str(threads),
                str(vcf),
            ]
        )
