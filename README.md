# SV-xQTL

## Structural variant–centric analysis of xQTLs across molecular phenotypes

### Overview
1. Create SV genotyping reference panel
2. Run idxdepth on sample crams and parse json to sample manifest
3. Run genotyping with paragraph
4. Filter genotypes for QC, HWE, and call fraction

## Installation
This package is currently in active developlment. Install it in editable mode via:
```sh
git clone https://github.com/sciencesteveho/sveqtl
cd sveqtl
pip install -e .
```

Additionally, bcftools and htslib installations in the environment PATH are required.


## Examples
Build the SV reference for genotyping:
```shell
python -u "${genotype_script}" \
    --output_path "${outdir}" \
    --reference "${reference}"
```