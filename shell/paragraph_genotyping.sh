#!/bin/bash

# TODO : Update M, threads, and initial vars. Fix any hardcoded dirs (like bin/)

sample=$1
sv_genotype_reference=/ocean/projects/bio210019p/stevesho/sveqtl/genotyping/sv_genotype_reference.vcf
manifest_dir=/ocean/projects/bio210019p/stevesho/sveqtl/genotyping/paragraph/genotype_resources
reference=/ocean/projects/bio210019p/stevesho/resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
geuvadis_pilot=/ocean/projects/bio210019p/stevesho/sveqtl/genotyping/geuvadis_pilot


# =============================================================================
# Filter genotype calls removing errant terms from genotyping QC
# Args:
#   1 - Directory containing the pilot data
#   2 - Sample name
# Returns:
#   Writes filtered genotypes to a gzipped VCF file
# Examples::
#   filter_genotypes /path/to/pilot_data sample_name
# =============================================================================
filter_genotypes() {
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
# Run Paragraph genotyping for a given sample
# Args:
#   1 - Sample name
# Returns:
#   Writes raw genotypes.vcf.gz into the sample directory
# =============================================================================
paragraph_genotyping() {
  local sample=$1
  local out_dir="${geuvadis_pilot}/${sample}"
  mkdir -p "$out_dir"

  python -u bin/multigrmpy.py \
    -i "$sv_genotype_reference" \
    -m "${manifest_dir}/${sample}.manifest.tsv" \
    -r "$reference" \
    -o "$out_dir" \
    -M 153 \
    -t 32
}

# =============================================================================
# Index the filtered VCF with tabix
# Args:
#   1 - Directory containing the pilot data
#   2 - Sample name
# Returns:
#   Creates a .tbi index for the filtered VCF
# =============================================================================
index_vcf() {
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
  paragraph_genotyping "$sample"

  # filter genotypes
  filter_genotypes "$geuvadis_pilot" "$sample"

  # index the filtered VCF
  index_vcf "$geuvadis_pilot" "$sample"
}

# Run the main function
main "$sample"