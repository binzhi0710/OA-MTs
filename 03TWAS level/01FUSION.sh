#!/bin/bash

# Set directories
sumstats_dir="./MTAGGWAS"
weights_dir="./GTExWEIGHTS"
ld_ref_prefix="./LDREFEUR/1000G.EUR."
output_dir="./fusion_results"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Get all .sumstats and .pos files
sumstats_files=("$sumstats_dir"/*.sumstats)
weights_files=("$weights_dir"/*.pos)

# Iterate over chromosomes (1 to 22)
for chr in {1..22}
do
  # Iterate over sumstats files
  for sumstats_file in "${sumstats_files[@]}"
  do
    sumstats_base=$(basename "$sumstats_file" .sumstats)
    # Iterate over weights files
    for weights_file in "${weights_files[@]}"
    do
      weights_base=$(basename "$weights_file" .pos)
      # Define output filename
      output_file="$output_dir/${sumstats_base}_${weights_base}_chr${chr}.dat"
      
      # Generate and submit SLURM job
      sbatch <<EOF
#!/bin/bash
#SBATCH -J fusion_${sumstats_base}_${weights_base}_chr${chr}
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=15G
#SBATCH -p cyhpc_1
#SBATCH -t 20-12:00:00
#SBATCH --comment=test

Rscript FUSION.assoc_test.R \\
  --sumstats "$sumstats_file" \\
  --weights "$weights_file" \\
  --weights_dir "$weights_dir/" \\
  --ref_ld_chr "$ld_ref_prefix" \\
  --chr $chr \\
  --out "$output_file"
EOF
    done
  done
done

echo "All SLURM job scripts submitted."