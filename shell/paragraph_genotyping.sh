#!/bin/bash

"""Run Paragraph genotyping for a given sample and filter the genotype calls. Ensure environment has tabix available in the PATH."""


sample=$1

base_dir=/ocean/projects/bio210019p/stevesho/sveqtl/genotyping
sv_genotype_reference="${base_dir}"/sv_genotype_reference.vcf
manifest_dir="${base_dir}"/paragraph/genotype_resources
reference=/ocean/projects/bio210019p/stevesho/resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
geuvadis_pilot="${base_dir}"/geuvadis_pilot
idxdepth_dir="${base_dir}"/idxdepth
paragraph_dir="${base_dir}"/paragraph
threads=32


# =============================================================================
# Filter genotype calls removing errant terms from genotyping QC
# Args:
#   1 - Directory containing the pilot data
#   2 - Sample name
# Returns:
#   Writes filtered genotypes to a gzipped VCF file
# Examples::
#   _filter_genotypes /path/to/pilot_data sample_name
# =============================================================================
_filter_genotypes() {
  local pilot_dir=$1
  local sample=$2
  local bad_terms='BP_DEPTH|NO_READS|DEPTH|MULTIMATCHED|UNMATCHED|NO_VALID_GT|GQ|CONFLICT|BP_NO_GT'
  local infile="${pilot_dir}/${sample}/genotypes.vcf.gz"
  local outfile="${pilot_dir}/${sample}/genotypes_filtered.vcf.gz"

  zcat "$infile" \
    | awk -v re="$bad_terms" '
        /^#/ { print; next }
        !($0 ~ re)
      ' \
    | bgzip -c > "$outfile"
}


# =============================================================================
# Get M (max read alignments) from sample manifest
# Args:
#   1 - Full path to sample manifest file
# Returns:
#   20 * depth
# =============================================================================
_get_M() {
  local manifest=$1
  local depth

  depth=$(awk 'NR==2 {print int($3)}' "$manifest")
  echo $(( depth * 20 ))
}


# =============================================================================
# Run Paragraph genotyping for a given sample
# Args:
#   1 - Sample name
#   2 - Paragraph directory
#   3 - Manifest directory
#   4 - Out dir
#   5 - Threads to use (default: 32)
# Returns:
#   Writes raw genotypes.vcf.gz into the sample directory
# =============================================================================
_paragraph_genotyping() {
  local sample=$1
  local paragraph_dir=$2
  local manifest_dir=$3
  local geuvadis_pilot=$4
  local threads=$5

  local out_dir="${geuvadis_pilot}/${sample}"
  local manifest="${manifest_dir}/${sample}.manifest.tsv"

  local M
  mkdir -p "$out_dir"

  M=$(_get_M "$manifest")

  python -u "${paragraph_dir}"/bin/multigrmpy.py \
    -i "$sv_genotype_reference" \
    -m "${manifest_dir}/${sample}.manifest.tsv" \
    -r "$reference" \
    -o "$out_dir" \
    -M "$M" \
    -t "$threads"
}


# =============================================================================
# Index the filtered VCF with tabix
# Args:
#   1 - Directory containing the pilot data
#   2 - Sample name
# Returns:
#   Creates a .tbi index for the filtered VCF
# =============================================================================
_index_vcf() {
  local pilot_dir=$1
  local sample=$2

  tabix -p vcf "${pilot_dir}/${sample}/genotypes_filtered.vcf.gz"
}


# =============================================================================
# Main function to perform centralized processing
# =============================================================================
main() {
  local sample=$1

  # make required directory
  mkdir -p "${geuvadis_pilot}"/"${sample}"

  # genotype with paragraph
  _paragraph_genotyping \
    "$sample" \
    "$paragraph_dir" \
    "$manifest_dir" \
    "$geuvadis_pilot" \
    "$threads"

  # filter genotypes
  _filter_genotypes \
    "$geuvadis_pilot" \
    "$sample"

  # index the filtered VCF
  _index_vcf \
    "$geuvadis_pilot" \
    "$sample"
}


# Run the main function
main "$sample"