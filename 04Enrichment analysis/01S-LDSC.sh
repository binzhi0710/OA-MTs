#!/bin/bash

# Set paths
input_dir="./metaOAprepared"
output_dir="./metaOAresultV1.2"
ref_ld_chr="1000G_Phase3_baseline_v1.2_ldscores/baseline_v1.2/baseline."
frqfile_chr="1000G_Phase3_frq/1000G.EUR.QC."
w_ld_chr="1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."

# Create output directory
mkdir -p "$output_dir"

# Set maximum parallel jobs
max_jobs=10
running_jobs=0

# Process each .sumstats.gz file
for sumstats_file in "$input_dir"/*.sumstats.gz; do
    filename=$(basename "$sumstats_file" .sumstats.gz)
    result_file="$output_dir/${filename}.log"

    # Skip completed files
    if [[ -f "$result_file" ]]; then
        echo "Skipping completed file: $filename"
        continue
    fi

    # Run LDSC analysis in background
    /usr/bin/python2.7 ldsc.py \
        --h2 "$sumstats_file" \
        --ref-ld-chr "$ref_ld_chr" \
        --frqfile-chr "$frqfile_chr" \
        --w-ld-chr "$w_ld_chr" \
        --overlap-annot \
        --print-cov \
        --print-coefficients \
        --print-delete-vals \
        --out "$output_dir/${filename}" &

    # Manage parallel jobs
    running_jobs=$((running_jobs + 1))
    if [[ $running_jobs -ge $max_jobs ]]; then
        wait
        running_jobs=0
    fi
done

# Wait for all jobs to complete
wait
echo "All tasks completed"