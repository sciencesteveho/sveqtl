#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""Merge SV callsets to produce a genotype reference panel."""


from typing import Dict, List

from sveqtl.catalogue_config import SOURCE_CONFIGS
from sveqtl.ref_downloader import RefDownloader
from sveqtl.sv_catalogue import SVCatalogue
from sveqtl.sv_reference_merger import SVReferenceMerger


def main() -> None:
    """Produce an SV genotyping reference panel."""
    short_read_svs: List[Dict] = []
    long_read_svs: List[Dict] = []

    # Download and preprocess SV reference callsets
    downloader = RefDownloader(configs=SOURCE_CONFIGS)
    downloader.download_ref_callsets()
    downloader.filter_callsets()

    # Normalize each VCF
    for source_name, config in SOURCE_CONFIGS.items():
        filtered_vcf_path = f"{config.output_name}.filtered.vcf"

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
    merger.write_merged("sv_genotype_reference.vcf")
    print("[RefMerger] Finished merging SV callsets.")


if __name__ == "__main__":
    main()
