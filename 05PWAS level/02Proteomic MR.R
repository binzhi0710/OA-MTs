library(TwoSampleMR)
library(ggplot2)
library(foreach)
library(data.table)

# Set directories and exposure files
expo_files <- c("decode.txt", "ukbppp.txt", "ARIC.txt")
outcome_files <- list.files(path = "./1030MTAGresult", pattern = "\\.txt$", full.names = TRUE)

# Create output directory
dir.create("1030MRresult", showWarnings = FALSE)

# Process each exposure file
for (expo_file in expo_files) {
  expo_rt <- fread(expo_file)
  expo_name <- tools::file_path_sans_ext(basename(expo_file))
  
  # Process each outcome file
  for (outcome_file in outcome_files) {
    outc_rt <- read_outcome_data(
      snps = expo_rt$SNP,
      filename = outcome_file,
      sep = "\t",
      snp_col = "SNP",
      beta_col = "mtag_beta",
      se_col = "mtag_se",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      eaf_col = "meta_freq",
      pval_col = "mtag_pval"
    )
    
    # Harmonize data
    harm_rt <- harmonise_data(
      exposure_dat = expo_rt, 
      outcome_dat = outc_rt,
      action = 2
    )
    harm_rt$id <- outcome_file
    
    # Perform MR analysis
    mr_result <- mr(harm_rt, method_list = c("mr_wald_ratio", "mr_ivw"))
    
    # Generate OR results
    result_or <- generate_odds_ratios(mr_result)
    
    # Save results and harmonized data
    if (!is.null(result_or)) {
      result_or$id <- outcome_file
      result_filename <- paste0("1030MRresult/", tools::file_path_sans_ext(basename(outcome_file)), "_", expo_name, "_mr_results.txt")
      write.table(result_or, result_filename, sep = "\t", quote = FALSE, row.names = FALSE)
      
      harmonized_filename <- paste0("1030MRresult/", tools::file_path_sans_ext(basename(outcome_file)), "_", expo_name, "_harmonized_data.txt")
      write.table(harm_rt, harmonized_filename, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
}
