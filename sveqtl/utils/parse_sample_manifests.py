#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""Create sample manifest for Paragraph from idxdepth JSON."""


import argparse
import glob
import json
import os
from typing import Dict, List


def extract_sample_id_from_path(bam_path: str) -> str:
    """Extract sample ID from BAM path filename."""
    filename = os.path.basename(bam_path)
    return filename.split(".")[0]


def process_json(json_path: str) -> Dict[str, str]:
    """Process a single JSON file and extract required information."""
    try:
        with open(json_path, "r") as f:
            data = json.load(f)

        bam_path = data.get("bam_path", "")
        autosome_depth = data.get("autosome", {}).get("depth")
        read_length = data.get("read_length")
        sample_id = extract_sample_id_from_path(bam_path)

        return {
            "id": sample_id,
            "path": bam_path,
            "depth": autosome_depth,
            "read_length": read_length,
        }

    except Exception as e:
        raise ValueError(
            f"Failed to process JSON file {json_path}. "
            "Ensure it contains the required fields."
        ) from e


def write_sample_manifest(sample: Dict[str, str], output_file: str) -> None:
    """Write a single sample manifest to a TSV file."""
    with open(output_file, "w") as f:
        f.write("id\tpath\tdepth\tread length\n")
        f.write(
            f"{sample['id']}\t{sample['path']}\t{sample['depth']}\t{sample['read_length']}\n"
        )


def parse_sample_manifests(
    idxdepth_out: str,
    sample: str,
) -> None:
    """Create individual sample manifests from JSON files."""
    out_dir = f"{idxdepth_out}/manifests"
    os.makedirs(out_dir, exist_ok=True)
    json_file = f"{idxdepth_out}/{sample}.json"

    if not json_file:
        raise ValueError(
            "No JSON file found. "
            "Ensure idxdepth has been run and the JSON file exists."
        )

    if sample_data := process_json(json_file):
        manifest_filename = f"{sample_data['id']}.manifest.tsv"
        output_path = f"{out_dir}/{manifest_filename}"

        write_sample_manifest(sample_data, output_path)

    print(f"Manifest files written to: {out_dir}")
