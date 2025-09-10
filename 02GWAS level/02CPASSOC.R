# Load R packages and source script
library(dplyr)
library(stringr)
library(data.table)
library(future)
setDTthreads(32)
source("./CPASSOC_FunctionSet.R")
snplist <- fread(input = paste0("plink_pruning.prune.in"), header = FALSE, data.table = FALSE)

# Read trait pair file
trait_pairs <- fread("1030traitpair.txt")

# Define output paths
out_path <- "./data_cpassoc"
dir.create(out_path, showWarnings = FALSE)
out_res_path <- "./result_cpassoc"
dir.create(out_res_path, showWarnings = FALSE)

# Loop through each row for analysis
for (i in 1:nrow(trait_pairs)) {
  
  # Read trait1 and trait2 information
  trait1 <- sprintf("%02d", as.numeric(trait_pairs$trait1[i]))
  trait2 <- trait_pairs$trait2[i]
  
  trait1_path <- file.path("./00metabo/dataprepared0.01maf/", paste0(trait1, ".txt"))
  trait2_path <- file.path("./00OA_noMHC/maf0.01/", paste0(trait2, ".txt"))
  
  # Read and format trait1 file
  data1 <- fread(trait1_path) %>%
    rename(SNP = SNP,
           effect_allele = A1,
           other_allele = A2,
           beta = beta,
           se = se,
           eaf = eaf,
           pval = pval,
           samplesize = N,
           z = Z) %>%
    mutate(chr = as.numeric(chr),
           pos = as.numeric(pos)) %>%
    na.omit()  # Remove rows with NA values
  
  # Read and format trait2 file
  data2 <- fread(trait2_path) %>%
    rename(SNP = SNP,
           effect_allele = A1,
           other_allele = A2,
           beta = beta,
           se = se,
           eaf = eaf,
           pval = pval,
           samplesize = N,
           z = Z) %>%
    mutate(chr = as.numeric(chr),
           pos = as.numeric(pos),
           effect_allele = toupper(effect_allele),  # Convert effect_allele to uppercase
           other_allele = toupper(other_allele)) %>%
    na.omit()  # Remove rows with NA values
  
  # Define paths for formatted trait1 and trait2 files
  formatted_trait1_path <- paste0(out_path, "/", trait1, "_formatted.txt")
  formatted_trait2_path <- paste0(out_path, "/", trait2, "_formatted.txt")
  
  # Save formatted data
  fwrite(data1, file = formatted_trait1_path, sep = "\t", quote = FALSE, row.names = FALSE)
  fwrite(data2, file = formatted_trait2_path, sep = "\t", quote = FALSE, row.names = FALSE)
  
  sample_size <- c(max(data1$samplesize), max(data2$samplesize))
  
  # Combine GWAS data
  d <- cpassoc.combine_assoc(files_path = c(formatted_trait1_path, formatted_trait2_path),
                             traits = c(trait1, trait2))
  
  # Define output file paths
  fn_out <- paste0(out_path, "/", trait1, "_", trait2, ".txt")
  fn_out_CPASSOC <- paste0(out_res_path, "/", trait1, "_", trait2, "_CPASSOC_SHet_SHom.txt")
  
  # Save combined data
  fwrite(d, file = fn_out, quote = FALSE, sep = "\t", row.names = FALSE)
  
  # Filter and calculate correlation matrix
  cor_mat <- d %>%
    # Filter 1: Select independent SNPs
    dplyr::filter(SNP %in% snplist$V1) %>%
    # Filter 2: Remove SNPs with z-score > 1.96
    dplyr::filter(dplyr::if_all(tidyselect::starts_with('z'), ~ abs(.x) < 1.96)) %>%
    dplyr::select(tidyselect::starts_with('z')) %>%
    as.matrix() %>%
    cor()
  
  # Perform SHet and SHom tests
  options(future.globals.maxSize = 4 * 1024^3)  # Set to 4GB
  
  d_output <- cpassoc.shet_shom_test(
    gwas_data = d,
    SHet_test = TRUE,
    SHom_test = TRUE,
    workers = 24,  # Adjust number of workers based on system resources
    mvrnorm_n = 1e6,
    set_seed = 666
  )
  
  # Save results
  fwrite(d_output, file = fn_out_CPASSOC, quote = FALSE, sep = "\t", row.names = FALSE)
  
  message("Completed analysis for: ", trait1, " vs ", trait2)
}

message("All analyses completed.")