library(readxl)
library(data.table)
library(hyprcoloc)

# Read Excel file
merged_data <- read_xlsx("merged_coloc_regions02.xlsx")

# Total number of loci
total_loci <- nrow(merged_data)

# Initialize list to store results
results_list <- list()

# Iterate through each row
for (j in 1:total_loci) {
  chr_locus <- merged_data$chr[j]
  start <- merged_data$start[j]
  end <- merged_data$end[j]
  
  # Extract non-NA file columns for the current row
  gwas_files <- merged_data[j, -c(1:3)] %>% unlist() %>% na.omit()
  
  # Initialize list to store filtered results for each GWAS file
  gwas_subsets <- list()
  snp_set <- NULL  # To store common SNPs across all files
  
  # Iterate through each GWAS file
  for (k in seq_along(gwas_files)) {
    file_path <- gwas_files[k]
    
    # Select only SNP, chr, pos, beta, and se columns
    gwas_data <- fread(file_path, select = c("SNP", "chr", "pos", "beta", "se"))
    
    # Filter SNPs for the specified chromosome and range
    gwas_subset <- gwas_data[chr == chr_locus & pos >= start & pos <= end, .(SNP, beta, se)]
    
    # Remove rows with se equal to 0
    gwas_subset <- gwas_subset[se != 0]
    
    # Ensure SNP column is unique, removing duplicate SNPs
    gwas_subset <- unique(gwas_subset, by = "SNP")
    
    # Store filtered results in the list
    gwas_subsets[[k]] <- gwas_subset
    
    # Update common SNP set
    if (is.null(snp_set)) {
      snp_set <- gwas_subset$SNP
    } else {
      snp_set <- intersect(snp_set, gwas_subset$SNP)
    }
  }
  
  # Filter common SNPs in each gwas_subset
  beta_list <- list()
  se_list <- list()
  
  for (k in seq_along(gwas_subsets)) {
    gwas_filtered <- gwas_subsets[[k]][SNP %in% snp_set]
    beta_list[[paste0("beta_", k)]] <- gwas_filtered$beta
    se_list[[paste0("se_", k)]] <- gwas_filtered$se
  }
  
  # Combine beta and se lists
  betas <- do.call(cbind, beta_list)
  ses <- do.call(cbind, se_list)
  
  # Set row names to SNPs
  rownames(betas) <- snp_set
  rownames(ses) <- snp_set
  
  # Perform hyprcoloc analysis
  traits <- paste0("T", seq_along(gwas_files))  # Set trait names
  res <- hyprcoloc(betas, ses, trait.names = traits, snp.id = snp_set)
  
  # Save results and add locus information
  res_df <- data.frame(do.call("cbind", res))
  res_df$chr <- chr_locus
  res_df$start <- start
  res_df$end <- end
  res_df$files <- paste(gwas_files, collapse = "; ")  # Combine file paths into a single string
  
  # Store results in the list
  results_list[[paste0("locus_", j)]] <- res_df
  
  # Print progress message
  cat("Completed hyprcoloc analysis for locus", j, "of", total_loci, "loci,", total_loci - j, "loci remaining\n")
}

# Combine all results
combined_df <- do.call(rbind, results_list)

# Save results to CSV file
write.csv(combined_df, "hyprcoloc_results.csv", row.names = TRUE)