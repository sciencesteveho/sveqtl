#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF >&2
Usage:
  $0 PARAGRAPH_DIR BAM_OR_CRAM_FILE REFERENCE_FASTA THREADS OUTPUT_DIR

Arguments:
  PARAGRAPH_DIR    Path to Paragraph install.
  BAM_OR_CRAM_FILE Input BAM/CRAM file to analyze.
  REFERENCE_FASTA  Reference FASTA.
  THREADS          Number of threads to use.
  OUTPUT_DIR       Directory where JSON output should be written.

Returns:
  The sample prefix (basename without extension) on success.
EOF
  exit 1
}

#######################################
# Run idxdepth for a single sample.
# Globals:
#   None
# Arguments:
#   paragraph_dir: Path to Paragraph installation directory
#   bam_file     : Input BAM/CRAM file path
#   reference    : Reference FASTA path
#   threads      : Number of threads to use
#   output_dir   : Directory where JSON output should be written
# Returns:
#   Prefix (basename without extension)
#######################################
run_idxdepth() {
  if [[ $# -lt 5 ]]; then
    usage
  fi

  local paragraph_dir="$1"
  local bam_file="$2"
  local reference="$3"
  local threads="$4"
  local output_dir="$5"
  local prefix

  prefix="$(basename "$bam_file" | cut -d. -f1)"

  mkdir -p "$output_dir"

  "${paragraph_dir}/bin/idxdepth" \
    --bam "$bam_file" \
    --reference "$reference" \
    --threads "$threads" \
    --output "${output_dir}/${prefix}.json"

  echo "$prefix"
}

# Invoke main
run_idxdepth "$@"
