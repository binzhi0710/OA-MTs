setwd("~/S-CARMA/")
library(data.table)
library(magrittr)
library(dplyr)
library(devtools)
library(R.utils)
library(arrow)
library(CARMA)

# Get list of files in 03plink directory
file_list <- list.files(path = "./03plink", pattern = "_plink.txt$", full.names = TRUE)

# Get processed files from CARMAresult directory
processed_files <- list.files(path = "./CARMAresult", pattern = ".txt.gz$", full.names = FALSE)
processed_base_names <- gsub(".txt.gz$", "", processed_files)

# Filter unprocessed files
file_list_to_process <- file_list[!basename(file_list) %in% paste0(processed_base_names, "_plink.txt")]

# Set number of files per script
files_per_script <- 1
total_files <- length(file_list_to_process)
total_scripts <- ceiling(total_files / files_per_script)

# Generate R scripts
for (i in 1:total_scripts) {
  script_name <- paste0("scriptOAMTs_", sprintf("%03d", i), ".R")
  
  # Determine file range for this script
  start_idx <- (i - 1) * files_per_script + 1
  end_idx <- min(i * files_per_script, total_files)
  current_files <- file_list_to_process[start_idx:end_idx]
  
  # Initialize script content
  script_content <- paste0(
    "setwd('~/S-CARMA/')\n",
    "library(data.table)\n",
    "library(magrittr)\n",
    "library(dplyr)\n",
    "library(devtools)\n",
    "library(R.utils)\n",
    "library(arrow)\n",
    "library(CARMA)\n",
    "\nplink_path <- path.expand('~/S-CARMA/plink/plink')\n",
    "bfile_path <- './EUR2m'\n\n"
  )
  
  # Generate processing code for each file
  for (snp_file in current_files) {
    base_name <- gsub("_plink.txt$", "", basename(snp_file))
    mtag_base_name <- gsub("_locus[0-9]+", "", base_name)
    
    file_processing_code <- paste0(
      "output_prefix <- './05plinked/", base_name, "_plink'\n",
      "ld_file <- './05plinked/", base_name, "_plink.ld'\n",
      "locus1polyfuned_file <- './04polyfuned/", base_name, "_polyfuned'\n",
      "parquet_data_file <- './02annotation/", base_name, "_annot.txt'\n",
      "mtag_file <- '~/02SI/", mtag_base_name, ".txt'\n",
      "output_file <- 'CARMAresult/", base_name, ".txt.gz'\n",
      "\n",
      "# Read files\n",
      "locus1polyfuned <- read.table(locus1polyfuned_file, header = TRUE, sep = '\\t')\n",
      "parquet_data_filtered <- fread(parquet_data_file)\n",
      "\n",
      "# Merge and filter SNPs\n",
      "merged_data <- merge(parquet_data_filtered, locus1polyfuned[, c('SNP', 'SNPVAR')], by = 'SNP', all.x = TRUE)\n",
      "rearranged_data <- merged_data[, c('SNPVAR', setdiff(names(merged_data), 'SNPVAR')), with = FALSE]\n",
      "colnames(rearranged_data)[colnames(rearranged_data) == 'SNPVAR'] <- 'PolyFun'\n",
      "rearranged_data <- rearranged_data[order(rearranged_data$BP), ]\n",
      "rearranged_data <- rearranged_data[complete.cases(rearranged_data), ]\n",
      "rownames(rearranged_data) <- NULL\n",
      "\n",
      "# Read MTAG file\n",
      "data <- fread(mtag_file)\n",
      "snps_to_filter <- rearranged_data$SNP\n",
      "sumstat_filtered <- data[data$SNP %in% snps_to_filter, ]\n",
      "setnames(sumstat_filtered, old = c('SNP', 'chr', 'pos', 'A1', 'A2', 'eaf', 'beta', 'se', 'pval'),\n",
      "         new = c('SNP', 'CHR', 'BP', 'A1', 'A2', 'eaf', 'beta', 'se', 'pval'))\n",
      "sumstat_filtered <- sumstat_filtered[order(sumstat_filtered$BP), ]\n",
      "rownames(sumstat_filtered) <- NULL\n",
      "\n",
      "# Write PLINK input file\n",
      "fwrite(sumstat_filtered[, c('SNP', 'A1')], file = './03plink/", base_name, "_plink.txt', sep = '\\t', col.names = FALSE, row.names = FALSE, quote = FALSE)\n",
      "\n",
      "# Run PLINK to generate LD matrix\n",
      "system2(plink_path,\n",
      "        args = c('--bfile', bfile_path,\n",
      "                 '--a1-allele', './03plink/", base_name, "_plink.txt', '2', '1', \"'#'\",\n",
      "                 '--extract', './03plink/", base_name, "_plink.txt',\n",
      "                 '--out', output_prefix,\n",
      "                 '--r',\n",
      "                 '--matrix',\n",
      "                 '--memory', '20480'))\n",
      "\n",
      "# Read LD file\n",
      "ld <- fread(file = ld_file)\n",
      "ld[is.na(ld)] <- 0\n",
      "\n",
      "# Perform CARMA analysis\n",
      "z.list <- list(sumstat_filtered$beta / sumstat_filtered$se)\n",
      "ld.list <- list(as.matrix(ld))\n",
      "lambda.list <- list(1)\n",
      "annot.list <- list(as.matrix(cbind(1, rearranged_data %>% select(-(SNP:A2)))))\n",
      "\n",
      "CARMA.results_annot <- CARMA(z.list, ld.list, w.list = annot.list, lambda.list = lambda.list, input.alpha = 0, outlier.switch = TRUE)\n",
      "\n",
      "# Initialize CS column\n",
      "sumstat_finish <- sumstat_filtered %>% mutate(PIP = CARMA.results_annot[[1]]$PIPs, CS = 0)\n",
      "\n",
      "# Update CS values for credible sets\n",
      "credible_sets <- CARMA.results_annot[[1]]$`Credible set`[[2]]\n",
      "for (cs_index in seq_along(credible_sets)) {\n",
      "  credible_snps <- credible_sets[[cs_index]]\n",
      "  if (length(credible_snps) > 0) {\n",
      "    sumstat_finish$CS[credible_snps] <- cs_index\n",
      "  }\n",
      "}\n",
      "\n",
      "# Save final results\n",
      "fwrite(x = sumstat_finish, file = output_file, sep = '\\t', quote = FALSE, na = 'NA', row.names = FALSE, col.names = TRUE, compress = 'gzip')\n\n"
    )
    
    script_content <- paste0(script_content, file_processing_code)
  }
  
  # Write script to file
  write(script_content, file = script_name)
}

# To submit all generated scripts, run the following command in the terminal:
# for script in *.sh; do sbatch $script; done