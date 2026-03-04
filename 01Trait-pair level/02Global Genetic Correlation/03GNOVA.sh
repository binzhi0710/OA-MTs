#!/bin/bash
#SBATCH -J gnova_script
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=25G
#SBATCH -p cydell_1,cydell_2,cydell_3,cyhpc_1
#SBATCH -t 20-12:00:00
#SBATCH --comment=test

# Set input and output folder paths
MTs_folder="./1030remakemetapreparedGNOVA"
OA_folder="./OApreparedGNOVA"
output_folder="./1030remakeGNOVAresults"

# Create output folder
mkdir -p $output_folder

# Get list of .sumstats.gz files in MTs_folder
MT_files=("$MTs_folder"/*.sumstats.gz)
OA_files=("$OA_folder"/*.gz)

# Initialize counter
count=0

# Calculate total number of iterations
total_files=$(( ${#MT_files[@]} * ${#OA_files[@]} ))

# Start nested loop
for mt_file in "${MT_files[@]}"; do
  mt_base=$(basename "${mt_file%.*}")
  for OA_file in "${OA_files[@]}"; do
    oa_base=$(basename "${OA_file%.*}")

    # Set output file path
    output_file="$output_folder/${mt_base}-${oa_base}.results"

    # Print progress bar
    ((count++))
    progress=$((count * 100 / total_files))
    echo -ne "Progress: $progress%    \r"

    # Run GNOVA analysis
    /usr/bin/python2.7 ./gnova.py "$mt_file" "$OA_file" \
      --bfile ./bfiles/eur_chr@_SNPmaf5 \
      --out "$output_file"
  done
done

echo -e "\nAnalysis completed!"