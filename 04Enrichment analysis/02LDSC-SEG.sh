#!/bin/bash

# Set paths
export PATH_TO_LDSC=/usr/bin/python2.7
export LDSC_SCRIPT=./ldsc.py
export REF_LD_CHR=./s-LDSC/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.
export W_LD_CHR=./s-LDSC/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.

# Function for LDSC analysis with Franke
run_ldsc_franke() {
  oamt_file=$1
  base_name=$(basename "$oamt_file" .sumstats.gz)
  output_file="./SEGout/${base_name}_Franke"
  "$PATH_TO_LDSC" "$LDSC_SCRIPT" \
    --h2-cts "$oamt_file" \
    --ref-ld-chr "$REF_LD_CHR" \
    --w-ld-chr "$W_LD_CHR" \
    --ref-ld-chr-cts ./Franke.ldcts \
    --out "$output_file"
}

# Function for LDSC analysis with GTEx
run_ldsc_gtex() {
  oamt_file=$1
  base_name=$(basename "$oamt_file" .sumstats.gz)
  output_file="./SEGout/${base_name}_GTEx"
  "$PATH_TO_LDSC" "$LDSC_SCRIPT" \
    --h2-cts "$oamt_file" \
    --ref-ld-chr "$REF_LD_CHR" \
    --w-ld-chr "$W_LD_CHR" \
    --ref-ld-chr-cts ./GTEx.ldcts \
    --out "$output_file"
}

export -f run_ldsc_franke
export -f run_ldsc_gtex

# Run both analyses in parallel with up to 11 tasks
parallel -j 11 run_ldsc_franke ::: ./metaOAprepared/*.sumstats.gz
parallel -j 11 run_ldsc_gtex ::: ./metaOAprepared/*.sumstats.gz