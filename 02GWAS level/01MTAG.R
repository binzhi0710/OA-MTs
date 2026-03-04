# Load necessary library
library(data.table)

# Read trait pair file
trait_pairs <- fread("1030traitpair.txt")

# Base script content
base_script <- "#!/bin/bash
#SBATCH -J mtag_analysis
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=15G
#SBATCH -p cydell_1,cydell_2,cydell_3,cyhpc_1
#SBATCH -t 20-12:00:00
#SBATCH --comment=mtag_test

# Start MTAG analysis
/usr/bin/python2.7 mtag.py \\
"

# Iterate through trait pairs to generate scripts
for (i in seq_len(nrow(trait_pairs))) {
  trait1 <- sprintf("%02d", as.numeric(trait_pairs$trait1[i]))
  trait2 <- trait_pairs$trait2[i]
  
  # Generate MTAG analysis command
  script_content <- paste0(
    base_script,
    "  --sumstats ./1030meta_prepare/", trait1, ".txt,./OAprepare/", trait2, ".txt \\\n",
    "  --out ./1030MTAGresult/", trait1, "_", trait2, " \\\n",
    "  --stream_stdout \\\n",
    "  --perfect_gencov \\\n",
    "  --equal_h2 \\\n",
    "  --force\n"
  )
  
  # Save script as .sh file, named 001.sh, 002.sh, etc., directly to working directory
  script_name <- sprintf("%03d.sh", i)
  writeLines(script_content, script_name)
}

cat("All MTAG scripts have been saved to the working directory!\n")

# To submit all generated scripts, run the following command in the terminal:
# for script in *.sh; do sbatch $script; done