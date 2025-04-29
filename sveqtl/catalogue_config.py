#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""Configs to store SV catalogue-specific processing rules."""


from dataclasses import dataclass
from dataclasses import field
from typing import Dict, Optional


@dataclass
class SourceRulesConfig:
    """Hold source-specific rules for SV processing.

    Attributes:
      trust_end: Whether to trust the END field in the VCF.
      trust_svlen: Whether to trust the SVLEN field in the VCF.
      ins_seq_mode: How to handle insertion sequences.
      exceptions: Optional dictionary of exceptions for specific SV types.
    """

    trust_end: bool
    trust_svlen: bool
    ins_seq_mode: str
    exceptions: Optional[Dict[str, Dict[str, bool]]] = field(default_factory=dict)


@dataclass
class SourceConfig:
    """Hold configuration for each SV callset source.

    Attributes:
      download_url: URL to download the VCF file.
      filter_mode: Mode for filtering the VCF (e.g., "info", "id").
      priority: Priority of the source.
      read_type: Type of reads (e.g., "long_read", "short_read").
      rules: SourceRulesConfig object with rules for processing.
      output_name: Name of the output file.
      zipped: Whether the file is zipped or not.
    """

    download_url: str
    filter_mode: str
    priority: int
    read_type: str
    rules: SourceRulesConfig
    output_name: str
    zipped: bool = False


SOURCE_CONFIGS: Dict[str, SourceConfig] = {
    # No END, has SVLEN for DEL
    # Insertion seq is stored in ALT (if ALT != <INS>)
    "hgsvc3_insdel": SourceConfig(
        download_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/variants_GRCh38_sv_insdel_alt_HGSVC2024v1.0.vcf.gz",
        filter_mode="info",
        priority=1,
        read_type="long_read",
        rules=SourceRulesConfig(
            trust_end=False,
            trust_svlen=True,
            ins_seq_mode="alt",
        ),
        output_name="hgsvc3_insdel.vcf",
    ),
    # No END, but has positive SVLEN
    # No insertions
    "hgsvc3_inv": SourceConfig(
        download_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/variants_GRCh38_sv_inv_sym_HGSVC2024v1.0.vcf.gz",
        filter_mode="alt",
        priority=1,
        read_type="long_read",
        rules=SourceRulesConfig(
            trust_end=False,
            trust_svlen=True,
            ins_seq_mode="none",
        ),
        output_name="hgsvc3_inv.vcf",
    ),
    # No END, no SVLEN
    # Insertion seq is stored in ALT (if ALT != <INS>)
    "ont1k": SourceConfig(
        download_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/v1.1/giggles-genotyping/giggles-genotyping_final-vcf.unphased.vcf.gz",
        filter_mode="id",
        priority=2,
        read_type="long_read",
        rules=SourceRulesConfig(
            trust_end=False,
            trust_svlen=False,
            ins_seq_mode="alt",
        ),
        output_name="ont1k.vcf",
    ),
    # Has END, but for INS or INV needs to be recalculated
    # Insertion seq is stored in ALT
    "ont100": SourceConfig(
        download_url="https://s3.amazonaws.com/1000g-ont/Gustafson_etal_2024_preprint_SUPPLEMENTAL/20240423_jasmine_intrasample_noBND_custom_suppvec_alphanumeric_header_JASMINE.vcf.gz",
        filter_mode="info",
        priority=3,
        read_type="long_read",
        rules=SourceRulesConfig(
            trust_end=True,
            trust_svlen=True,
            ins_seq_mode="alt",
            exceptions={
                "INS": {"trust_end": False},
                "INV": {"trust_end": False, "trust_svlen": False},
            },
        ),
        output_name="ont100.vcf",
    ),
    # Has valid END for all
    # Insertion sequences are in INFO["INSEQ"]
    "1kg30x": SourceConfig(
        download_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz",
        filter_mode="alt",
        priority=4,
        read_type="short_read",
        rules=SourceRulesConfig(
            trust_end=True,
            trust_svlen=False,
            ins_seq_mode="info_inseq",
        ),
        output_name="1kg30x.vcf",
    ),
    # Has valid END for all
    # No insertions
    "ccdg": SourceConfig(
        download_url="https://github.com/hall-lab/sv_paper_042020/raw/refs/heads/master/Supplementary_File_1.zip",
        filter_mode="info",
        priority=5,
        read_type="short_read",
        rules=SourceRulesConfig(
            trust_end=True,
            trust_svlen=False,
            ins_seq_mode="none",
        ),
        zipped=True,
        output_name="ccdg.vcf",
    ),
}
