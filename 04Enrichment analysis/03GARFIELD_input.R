library(data.table)
library(tools)

# Set data.table threads
setDTthreads(threads = 0)

# Set paths
mtag_path <- "~/S-MTAG/MTAGresult/"
dataprepared_path <- "~/07metaboOA/dataprepared/"
noMHC_path <- "~/07metaboOA/00OA_noMHC/"
output_dir <- "./garfield-data/pval/"

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Get all MTAG files
files <- list.files(path = mtag_path, pattern = "*.txt", full.names = TRUE)

# Process each file
for (file_path in files) {
  mtag_data <- fread(file_path, select = c("SNP", "CHR", "BP", "mtag_pval"))
  mtag_filtered <- mtag_data[mtag_pval < 1e-5]
  
  file_prefix <- sub("_mtag_meta$", "", file_path_sans_ext(basename(file_path)))
  
  # Load prepared and noMHC files
  prepared_file <- file.path(dataprepared_path, paste0(sub("_.*$", "", file_prefix), ".txt"))
  prepared_data <- fread(prepared_file, select = c("SNP", "pval"))
  setnames(prepared_data, "pval", "pval_meta")
  
  noMHC_file <- file.path(noMHC_path, paste0(sub("^[0-9]+_", "", file_prefix), ".txt"))
  noMHC_data <- fread(noMHC_file, select = c("SNP", "pval"))
  setnames(noMHC_data, "pval", "pval_OA")
  
  # Merge and filter data
  merged_data <- merge(mtag_filtered, prepared_data, by = "SNP")
  merged_data <- merge(merged_data, noMHC_data, by = "SNP")
  final_data <- merged_data[mtag_pval < pval_meta & mtag_pval < pval_OA]
  
  # Create subfolder
  subfolder <- file.path(output_dir, file_prefix)
  dir.create(subfolder, showWarnings = FALSE)
  
  # Split by chromosome and save
  final_data_list <- split(final_data, final_data$CHR)
  for (chr in names(final_data_list)) {
    chr_data <- final_data_list[[chr]][, .(BP, mtag_pval)][order(BP)]
    fwrite(chr_data, file.path(subfolder, paste0("chr", chr)), col.names = FALSE, sep = "\t")
  }
}

# To execute this script, run in the terminal:
# Rscript process_garfield_data.R