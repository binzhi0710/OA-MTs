#!/bin/bash

OA_folder="./OA"
MTs_folder="./MTs"
output_folder="./MTs-OA"

mkdir -p $output_folder

OA_files=("$OA_folder"/*.gz)

MTs_files=("$MTs_folder"/*.gz)

total_files=$(( ${#OA_files[@]} * ${#MTs_files[@]} ))

count=0

for OA_file in "${OA_files[@]}"; do
  for MTs_file in "${MTs_files[@]}"; do
    OA_base=$(basename "${OA_file%.*}")
    MTs_base=$(basename "${MTs_file%.*}")

    output_file="$output_folder/${OA_base}-${MTs_base}.results"

    ((count++))
    progress=$((count * 100 / total_files))
    echo -ne "Progress: $progress%    \r"

    ./ldsc.py --rg "$OA_file,$MTs_file" --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out "$output_file"
  done
done

echo -e "\nAnalysis completed!"