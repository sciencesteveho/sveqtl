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


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract sample information from JSON files"
    )
    parser.add_argument(
        "--json_dir",
        help="Directory containing JSON files",
        default="/ocean/projects/bio210019p/stevesho/sveqtl/genotyping/idxdepth",
    )
    parser.add_argument(
        "--output_dir",
        default="/ocean/projects/bio210019p/stevesho/sveqtl/genotyping/paragraph/genotype_resources",
    )
    return parser.parse_args()


def _find_json_files(json_dir: str, pattern: str = "*.json") -> List[str]:
    """Find all JSON files in the specified directory matching the pattern."""
    json_pattern = os.path.join(json_dir, pattern)
    return glob.glob(json_pattern)


def write_sample_manifest(sample: Dict[str, str], output_file: str) -> None:
    """Write a single sample manifest to a TSV file."""
    with open(output_file, "w") as f:
        f.write("id\tpath\tdepth\tread length\n")
        f.write(
            f"{sample['id']}\t{sample['path']}\t{sample['depth']}\t{sample['read_length']}\n"
        )


def main() -> None:
    """Create individual sample manifests from JSON files."""
    args = parse_arguments()
    os.makedirs(args.output_dir, exist_ok=True)
    json_files = _find_json_files(args.json_dir)

    if not json_files:
        print("No JSON files found")
        return

    processed_count = 0

    for json_file in sorted(json_files):
        print(f"Processing: {os.path.basename(json_file)}")
        try:
            if sample_data := process_json(json_file):
                manifest_filename = f"{sample_data['id']}.manifest.tsv"
                output_path = f"{args.output_dir}/{manifest_filename}"

                write_sample_manifest(sample_data, output_path)
                processed_count += 1

        except Exception as e:
            raise ValueError(f"Error processing {json_file}: {e}") from e

    print(f"Successfully processed {processed_count} samples")
    print(f"Manifest files written to: {args.output_dir}")


if __name__ == "__main__":
    main()
