#!/bin/bash

#SBATCH -J PoPS
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --mem=60G
#SBATCH -p cyhpc_1,cydell_1,cydell_2,cydell_3
#SBATCH -t 20-12:00:00
#SBATCH --comment=test

out_files=$(ls *.genes.out)
module load apps/python/3.7.10

# Process each file
process_file() {
  file=$1
  prefix=$(echo $file | sed 's/\.genes.out$//')
  python3 pops.py \
    --gene_annot_path gene_all.txt \
    --feature_mat_prefix features_munged/pops_features \
    --num_feature_chunks 2 \
    --magma_prefix ./$prefix \
    --control_features_path features_jul17_control.txt \
    --out_prefix ./$prefix
}

export -f process_file

# Run concurrently with up to 5 processes
echo "$out_files" | xargs -n 1 -P 5 -I {} bash -c 'process_file "$@"' _ {}