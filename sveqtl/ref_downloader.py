#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""Download and preprocess SV reference callsets.

Note that this script is largely ad hoc due to specific differences between SV
VCFs, download URLS, and differing consortium conventions.
"""


from pathlib import Path
import subprocess
from typing import Dict

from sveqtl.catalogue_config import SourceConfig


class RefDownloader:
    """Class to download and prepare SV reference callsets.

    Our approach utilizes a mix of high-confidence, long-read, multi-platform
    resolved SVs, with population-scale short-read catalogues which help
    overcome the relatively small number of discovery sets from the
    multi-platform and improve variant diversity for genotyping. The reference
    callsets are as follows:

      1. Logsdon et al. - Multi-platform assemblies in 65 samples
      2. Schloissnig et al. - ONT in 1,019 1KG samples (callset is merged with
         HPRC1)
      3. Gustafson et al. - ONT in 100 1KG samples
      4. Byrska-Bishop et al. - Short read ~30X in 3,202 1KG samples
      5. Abel at al., Nature, 2020  - Short read ~20x in > 17,795 samples

    Processing steps:
     * Download reference callsets.
     * Remove sample information.
     * Filter callsets to keep PASS.
     * Filter the callsets to keep DEL, INS, DUP, and INV.
     * Collapse phase information, if present.

    Attributes:
        REF_CALLSETS: Dictionary mapping callset names to their download URLs.
        data_dir: Directory to store downloaded and processed files.

    Examples::
        # Instantiate the RefDownloader class and download reference callsets
        >>> downloader = RefDownloader(
        ...     output_path="path/to/output",
        ...     configs=SOURCE_CONFIGS
        ... )
        >>> downloader.download_ref_callsets()
        >>> downloader.filter_callsets()
    """

    def __init__(self, configs: Dict[str, SourceConfig], output_path: str = "") -> None:
        """Initialize the RefDownloader class."""
        self.configs = configs
        self.data_dir = Path(output_path)
        self.data_dir.mkdir(exist_ok=True)

    def download_ref_callsets(self) -> None:
        """Download reference callsets from the specified URLs."""
        for name, conf in self.configs.items():
            if conf.zipped:
                self._process_ccdg(name, conf)
            else:
                self._process_regular_callset(name, conf)

    def filter_callsets(self) -> None:
        """Clean up the callsets.

        Loop over every .vcf in self.data_dir and:
        - Remove sample-level info to collapse any phasing data
        - Keep only PASS records
        - Keep only SVTYPE in [DEL, INS, DUP, INV]

        The output is written to a new file <original>.filtered.vcf
        """
        for _, conf in self.configs.items():
            uncompressed_path = self.data_dir / conf.output_name
            if not uncompressed_path.exists():
                print(f"[Warning] {uncompressed_path} not found; skipping filter.")
                continue

            filtered_path = uncompressed_path.with_suffix(".vcf.filtered.vcf")
            self._filter_vcf(uncompressed_path, filtered_path, conf.filter_mode)

    def _process_regular_callset(self, name: str, cfg: SourceConfig) -> None:
        """Download a single .vcf.gz callset (if needed), decompress it, rename."""
        compressed_path = self.data_dir / f"{name}.vcf.gz"
        final_vcf = self.data_dir / cfg.output_name

        if final_vcf.exists():
            print(f"[Download] {final_vcf.name} already exists, skipping.")
            return

        if not compressed_path.exists():
            print(f"[Download] Downloading {name} -> {compressed_path.name}")
            self._download_file(cfg.download_url, str(compressed_path))

        if compressed_path.exists():
            self._decompress_gzip(str(compressed_path))

        old_vcf_path = self.data_dir / f"{name}.vcf"
        if old_vcf_path.exists() and old_vcf_path != final_vcf:
            old_vcf_path.rename(final_vcf)

    def _process_ccdg(self, name: str, cfg: SourceConfig) -> None:
        """All CCDG processing steps."""
        final_vcf = self.data_dir / cfg.output_name
        if final_vcf.exists():
            print(f"[Download] {final_vcf.name} already exists, skipping.")
            return

        zip_path = self.data_dir / f"{name}.zip"
        self._download_ccdg(cfg.download_url, zip_path)
        self._extract_ccdg(zip_path)
        self._process_ccdg_vcf(final_vcf)
        self._cleanup_ccdg_files(zip_path)

    def _download_ccdg(self, ccdg_url: str, zip_path: Path) -> None:
        """Download CCDG zip file."""
        if not zip_path.exists():
            self._download_file(ccdg_url, str(zip_path))

    def _extract_ccdg(self, zip_path: Path) -> None:
        """Extract contents from CCDG zip file."""
        subprocess.run(
            ["unzip", "-o", str(zip_path), "-d", str(self.data_dir)], check=True
        )

    def _process_ccdg_vcf(self, final_vcf_path: Path) -> None:
        """Process the extracted CCDG VCF file."""
        extracted_vcf_gz = self.data_dir / "Build38.public.v2.vcf.gz"
        if extracted_vcf_gz.exists():
            self._decompress_gzip(str(extracted_vcf_gz))
            decompressed_path = self.data_dir / "Build38.public.v2.vcf"
            if decompressed_path.exists():
                decompressed_path.rename(final_vcf_path)

    def _cleanup_ccdg_files(self, zip_path: Path) -> None:
        """Remove unnecessary CCDG files and rename the VCF."""
        bedpe_path = self.data_dir / "Build38.public.v2.bedpe.gz"
        if bedpe_path.exists():
            bedpe_path.unlink()
        if zip_path.exists():
            zip_path.unlink()

    def _filter_vcf(self, input_vcf: Path, output_vcf: Path, filter_mode: str) -> None:
        """Filter VCF files to keep only relevant variants."""
        base_filter = '(FILTER="PASS" || FILTER=".")'
        sv_filter = self._filter_vcf_by_mode(filter_mode)

        filter_expr = f"{base_filter} && {sv_filter}"
        bcftools_cmd = [
            "bcftools",
            "view",
            "-i",
            filter_expr,
            "-G",
            "-O",
            "v",
            "-o",
            str(output_vcf),
            str(input_vcf),
        ]
        subprocess.run(bcftools_cmd, check=True)

    @staticmethod
    def _filter_vcf_by_mode(filter_mode: str) -> str:
        """Return bcftools expression to keep specific SV types."""
        if filter_mode == "alt":
            return '(ALT="<DEL>" || ALT="<INS>" || ALT="<DUP>" || ALT="<INV>")'
        elif filter_mode == "id":
            return '(ID~"-DEL-" || ID~"-INS-" || ID~"-DUP-" || ID~"-INV-")'
        elif filter_mode == "info":
            return '(INFO/SVTYPE="DEL" || INFO/SVTYPE="INS" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV")'
        else:
            raise ValueError(f"Unknown filter_mode: {filter_mode}")

    @staticmethod
    def _download_file(url: str, output_path: str) -> None:
        """Download a file from URL to output path."""
        subprocess.run(["wget", "-O", output_path, url], check=True)

    @staticmethod
    def _decompress_gzip(file_path: str) -> None:
        """Decompress a gzipped file."""
        subprocess.run(["gunzip", file_path], check=True)
