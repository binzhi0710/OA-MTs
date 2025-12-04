#!/bin/bash
#SBATCH -J task_x9k
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=25G
#SBATCH -p node_a,node_b,node_c,node_d
#SBATCH -t 20-12:00:00
#SBATCH --comment=run

MTs_folder="./1030remakemetapreparedGNOVA"
OA_folder="./OApreparedGNOVA"
output_folder="./1030remakeGNOVAresults"

mkdir -p $output_folder

MT_files=("$MTs_folder"/*.sumstats.gz)
OA_files=("$OA_folder"/*.gz)

for mt_file in "${MT_files[@]}"; do
  mt_base=$(basename "${mt_file%.*}")
  for OA_file in "${OA_files[@]}"; do
    oa_base=$(basename "${OA_file%.*}")
    output_file="$output_folder/${mt_base}-${oa_base}.results"
    /usr/bin/python2.7 ./gnova.py "$mt_file" "$OA_file" \
      --bfile ./bfiles/eur_chr@_SNPmaf5 \
      --out "$output_file"
  done
done
