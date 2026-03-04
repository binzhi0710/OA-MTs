# Set folder paths
gwas_folder <- "./MTAG"
output_folder <- "./"
result_folder <- "./mbatresult/"

# Create result folder if it doesn't exist
if (!dir.exists(result_folder)) {
  dir.create(result_folder)
}

# Get all .ma files
gwas_files <- list.files(gwas_folder, pattern = "\\.ma$", full.names = TRUE)

# Base script content
base_script <- "#!/bin/bash
#SBATCH -J gcta_analysis
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=15G
#SBATCH -p cyhpc_1
#SBATCH -t 20-12:00:00
#SBATCH --comment=test

# Add executable permission to gcta64
chmod +x gcta64

"

# Generate shell scripts
for (i in seq_along(gwas_files)) {
  gwas_file <- gwas_files[i]
  gwas_base <- tools::file_path_sans_ext(basename(gwas_file))
  
  # Generate script content for each file
  script_content <- paste0(
    base_script,
    "./gcta64 --bfile g1000eur \\\n",
    "         --mBAT-combo ", gwas_file, " \\\n",
    "         --mBAT-gene-list glist_ensgid_hg19_v40.txt \\\n",
    "         --out ", result_folder, gwas_base, " \\\n",
    "         --diff-freq 0.1 \\\n",
    "         --thread-num 10\n"
  )
  
  # Save script as .sh file (e.g., 001.sh, 002.sh)
  script_name <- sprintf("%03d.sh", i)
  writeLines(script_content, file.path(output_folder, script_name))
}

cat("All scripts generated successfully!\n")

# To execute all generated scripts, run the following command in the terminal:
# for script in *.sh; do sbatch $script; done