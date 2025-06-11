# SV-xQTL

## Structural variant–centric analysis of xQTLs across molecular phenotypes

### Overview
1. Create SV genotyping reference panel
2. Run idxdepth on sample crams and parse json to sample manifest
3. Run genotyping with paragraph
4. Filter genotypes for QC, HWE, and call fraction

</br>

## Installation
This package is currently in active developlment. Install it in editable mode via:
```sh
git clone https://github.com/sciencesteveho/sveqtl
cd sveqtl
pip install -e .
```

Additionally, ensure that bcftools and htslib are installed in your environment’s PATH.

</br>

## Usage


### 1. Build sv genotype reference

```bash
sveqtl build_genotype_ref \
  --reference /path/to/reference.fasta \
  --output_dir /path/to/output/
```

**Required arguments:**

* `--reference`: Path to reference FASTA file (index must be in same directory)
* `--output_dir`: Output location for the SV genotyping reference panel

**Optional arguments:**

* `--concordance`: Only keep variants supported by 2+ datasets

### 2. Run idxdepth genotyping prerequisite

```bash
sveqtl run_idx_depth \
  --paragraph_dir /path/to/paragraph/ \
  --bam_file /path/to/sample.bam \
  --reference /path/to/reference.fasta \
  --threads 8 \
  --idxdepth_out /path/to/output/
```

**Required arguments:**

* `--paragraph_dir`: Path to Paragraph installation
* `--bam_file`: Input CRAM/BAM file to process
* `--reference`: Reference FASTA file
* `--idxdepth_out`: Output directory for results

**Optional arguments:**

* `--threads`: Number of threads to use (default: 4)

### 3. Genotype SVs using paragraph

```bash
sveqtl genotype \
  --paragraph_dir /path/to/paragraph/ \
  --manifest_dir /path/to/manifests/ \
  --reference /path/to/reference.fasta \
  --pilot_dir /path/to/pilot_output/ \
  --threads 32 \
  SAMPLE_NAME
```

**Required arguments:**

* `--paragraph_dir`: Path to Paragraph installation directory
* `--manifest_dir`: Directory containing per-sample manifest TSVs
* `--reference`: Reference FASTA file
* `--pilot_dir`: Output directory for genotyping results
* `SAMPLE_NAME`: Sample name (must match manifest basename)

**Optional arguments:**

* `--threads`: Number of threads to use (default: 32)

### 4. Combine per-sample VCFs into cohort and apply HWE + call fraction filters

```bash
sveqtl combine_and_filter_genotypes \
  --sample_dir /path/to/sample_vcfs/ \
  --threads 18 \
  --max_missing_frac 0.5 \
  --hwe_alpha 1e-4 \
  --min_maf 0.05
```

**Required arguments:**

* `--sample_dir`: Directory containing per-sample VCF files

**Optional arguments:**

* `--threads`: Number of threads to use (default: 18)
* `--max_missing_frac`: Maximum allowed missing genotype fraction (default: 0.50)
* `--hwe_alpha`: Hardy-Weinberg equilibrium significance threshold (default: 1e-4)
* `--min_maf`: Minimum minor allele frequency (default: 0.05)

  \--hwe\_alpha 1e-4&#x20;
  \--min\_maf 0.05
