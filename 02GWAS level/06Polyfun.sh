#!/bin/bash

# Input and output directories
input_dir="./01polyfun"
output_dir="./04polyfuned"

# Create output directory if it doesn't exist
mkdir -p $output_dir

# Get all .txt files in input directory
files=($input_dir/*.txt)

# Calculate total number of files
total_files=${#files[@]}

# Process 15 files per job
batch_size=15

# Submit SLURM jobs for each batch
for ((i=0; i<$total_files; i+=$batch_size))
do
  batch_files=("${files[@]:$i:$batch_size}")
  file_names=$(printf "%s " "${batch_files[@]}")
  
  sbatch <<EOF
#!/bin/bash
#SBATCH -J polyfun_batch_$((i / batch_size + 1))
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=8G
#SBATCH -p cyhpc_1
#SBATCH -t 20-12:00:00
#SBATCH --comment=test

module load apps/python/3.7.10

# Process each file in the batch
for input_file in $file_names
do
  base_name=\$(basename \$input_file .txt)
  output_file="$output_dir/\${base_name}ed"
  python3 ./polyfun/extract_snpvar.py --sumstats \$input_file --out \$output_file
  echo "Processed \$input_file, output to \$output_file"
done
EOF
done