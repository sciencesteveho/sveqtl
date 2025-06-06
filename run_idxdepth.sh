#!/bin/bash

# =============================================================================
# Run idxdepth for a given BAM/CRAM file
# Args:
#   1 - BAM/CRAM file
#   2 - Reference genome (e.g., GRCh38)
#   3 - Number of threads (default: 16)
#   4 - Output directory (default: current directory)
# Returns:
#   0 on success
#   1 on failure
# Examples::
#    run_idxdepth \
#        /ocean/projects/bio210019p/stevesho/sveqtl/genotyping/1kgp/HG00733.final.cram \
#        /ocean/projects/bio210019p/stevesho/resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
#        16 \
#        /ocean/projects/bio210019p/stevesho/sveqtl/genotyping/1kgp
# =============================================================================
run_idxdepth() {
    local bam_file=$1
    local reference=$2
    local threads=${3:-16}
    local output_dir=${4:-./}
    local prefix=$(basename "$bam_file" | cut -d. -f1)

    mkdir -p "$output_dir"

    /ocean/projects/bio210019p/stevesho/sveqtl/genotyping/paragraph/bin/idxdepth \
        --bam "$bam_file" \
        --reference "$reference" \
        --threads "$threads" \
        --output "$output_dir/${prefix}.json"
}