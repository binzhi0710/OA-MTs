library(stringr)

# Set folder paths
gwas_folder <- "./MTAGGWAS/"
output_folder <- "./scripts/"
gtexv8_folder <- "./GTExv8/"

# Create output scripts directory
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# Get all .txt and .besd files
gwas_files <- list.files(gwas_folder, pattern = "\\.txt$", full.names = TRUE)
gtexv8_files <- list.files(gtexv8_folder, pattern = "\\.besd$", full.names = TRUE)

# Base script content
base_script <- "#!/bin/bash
#SBATCH -J rscript
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem=15G
#SBATCH -p cydell_1,cydell_2,cydell_3,cyhpc_1
#SBATCH -t 20-12:00:00
#SBATCH --comment=test

chmod +x smr

# eQTLGEN analysis
#./smr --bfile g1000eur --gwas-summary {gwas_file} --beqtl-summary ./00eQTLGEN/eQTL --out results/{gwas_base}_eQTLGEN --thread-num 10 --diff-freq-prop 0.95

# GTExv8 analysis
"

# Generate SMR commands for GTExv8 files
gtexv8_analysis <- function(gwas_file, gwas_base) {
  analysis <- ""
  for (gtexv8_file in gtexv8_files) {
    gtexv8_base <- tools::file_path_sans_ext(basename(gtexv8_file))
    analysis <- paste0(analysis, "./smr --bfile g1000eur --gwas-summary ", gwas_file, 
                       " --beqtl-summary ", gtexv8_file, 
                       " --out results/", gwas_base, "_", gtexv8_base, 
                       " --thread-num 10 --diff-freq-prop 0.95\n")
  }
  return(analysis)
}

# Generate shell scripts for each GWAS file
for (i in seq_along(gwas_files)) {
  gwas_file <- gwas_files[i]
  gwas_base <- tools::file_path_sans_ext(basename(gwas_file))
  
  # Replace placeholders in base script
  script_content <- gsub("{gwas_file}", gwas_file, base_script, fixed = TRUE)
  script_content <- gsub("{gwas_base}", gwas_base, script_content, fixed = TRUE)
  
  # Add GTExv8 analysis commands
  script_content <- paste0(script_content, gtexv8_analysis(gwas_file, gwas_base))
  
  # Save script as .sh file (e.g., 001.sh, 002.sh)
  script_name <- sprintf("%03d.sh", i)
  writeLines(script_content, file.path(script_name))
}

# Replace // with / in all .sh files
sh_files <- list.files(pattern = "\\.sh$")
for (file in sh_files) {
  file_content <- readLines(file)
  file_content <- gsub("//", "/", file_content)
  writeLines(file_content, file)
}

# Remove .besd from all .sh files
for (file in sh_files) {
  file_content <- readLines(file)
  file_content <- gsub("\\.besd", "", file_content)
  writeLines(file_content, file)
}

cat("All scripts generated and cleaned successfully!\n")

# To execute all generated scripts, run in the terminal:
# for script in *.sh; do sbatch "$script"; done