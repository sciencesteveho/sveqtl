#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF >&2
Usage: $0 SAMPLE_NAME SV_REF_VCF PARAGRAPH_DIR MANIFEST_DIR REFERENCE_FASTA GENOTYPE_OUT THREADS

  SAMPLE_NAME      Sample ID.
  SV_REF_VCF       Path to sv_genotype_reference.vcf.
  PARAGRAPH_DIR    Path to Paragraph install.
  MANIFEST_DIR     Directory of sample.manifest.tsv files.
  REFERENCE_FASTA  Reference FASTA.
  GENOTYPE_OUT     Output directory for genotypes.
  THREADS          Number of threads to use.
EOF
  exit 1
}

# Check args
if [[ $# -ne 7 ]]; then
  usage
fi

sample="$1"
sv_ref_vcf="$2"
paragraph_dir="$3"
manifest_dir="$4"
reference="$5"
genotype_out="$6"
threads="$7"

#######################################
# Compute M = 20 × (column 3 of line 2 in the manifest)
# Globals:
#   sample
#   manifest_dir
# Arguments:
#   None
# Outputs:
#   Echoes M (max read alignments) to stdout
#######################################
_get_M() {
  awk 'NR==2 { print int($3) * 20 }' "${manifest_dir}/${sample}.manifest.tsv"
}

#######################################
# Run Paragraph’s multigrmpy.py for a single sample
# Globals:
#   sample
#   sv_ref_vcf
#   paragraph_dir
#   manifest_dir
#   reference
#   genotype_out
#   threads
# Arguments:
#   None
# Outputs:
#   Writes raw genotypes VCF at "${genotype_out}/${sample}/genotypes.vcf.gz"
#######################################
_paragraph_genotyping() {
  local out_dir="${genotype_out}/${sample}"
  local manifest="${manifest_dir}/${sample}.manifest.tsv"
  local M

  mkdir -p "$out_dir"
  M="$(_get_M)"

  python -u "${paragraph_dir}/bin/multigrmpy.py" \
    -i "$sv_ref_vcf" \
    -m "$manifest" \
    -r "$reference" \
    -o "$out_dir" \
    -M "$M" \
    -t "$threads"
}

#######################################
# Remove uncallable genotypes from the VCF
# Globals:
#   sample
#   genotype_out
# Arguments:
#   None
# Outputs:
#   Writes filtered VCF at "${genotype_out}/${sample}/genotypes_filtered.vcf.gz"
#######################################
_remove_uncalled() {
  local in_vcf="${genotype_out}/${sample}/genotypes.vcf.gz"
  local out_vcf="${genotype_out}/${sample}/genotypes_filtered.vcf.gz"

  mkdir -p "${genotype_out}/${sample}"
  bcftools view -g ^miss -G -Oz -o "$out_vcf" "$in_vcf"
  bcftools index "$out_vcf"
}


#######################################
# Main entrypoint to drive genotyping pipeline
# Globals:
#   None
# Arguments:
#   None
#######################################
main() {
  _paragraph_genotyping
  _remove_uncalled
}

main "$@"
