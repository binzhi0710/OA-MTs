##### HDL Analysis #####
library(HDL)
library(data.table)
library(doSNOW)
setDTthreads(32)

# Set paths
LD.path <- "./LDimpute"
personality_file <- "./02HDLOA/AllOA.txt.gz"
pd_path <- "./00metabo/HDLpre"
output_dir <- "./02HDLresults"

# Get list of pd files
pd_files <- list.files(pd_path, pattern = "\\.txt.gz$", full.names = TRUE)

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Read the specified personality file
cat("Reading personality file: ", basename(personality_file), "\n")
gwas1 <- fread(personality_file)

# Loop through pd_files
for (d_file in pd_files) {
  cat("Processing d_file: ", basename(d_file), "\n")
  
  # Read d_file
  gwas2 <- fread(d_file)
  
  # HDL analysis
  res.HDL <- HDL.rg.parallel(gwas1, gwas2, LD.path, numCores = 4)
  
  # Create data frame
  res_df <- data.frame(
    rg = res.HDL$rg,
    rg_se = res.HDL$rg.se,
    P = res.HDL$P,
    Heritability_1 = res.HDL$estimates.df["Heritability_1", "Estimate"],
    Heritability_1_se = res.HDL$estimates.df["Heritability_1", "se"],
    Heritability_2 = res.HDL$estimates.df["Heritability_2", "Estimate"],
    Heritability_2_se = res.HDL$estimates.df["Heritability_2", "se"],
    Genetic_Covariance = res.HDL$estimates.df["Genetic_Covariance", "Estimate"],
    Genetic_Covariance_se = res.HDL$estimates.df["Genetic_Covariance", "se"],
    Genetic_Correlation = res.HDL$estimates.df["Genetic_Correlation", "Estimate"],
    Genetic_Correlation_se = res.HDL$estimates.df["Genetic_Correlation", "se"],
    eigen_use = res.HDL$eigen.use
  )
  
  # Build output file name
  p_file_name <- basename(sub("\\.txt.gz$", "", personality_file))
  d_file_name <- basename(sub("\\.txt.gz$", "", d_file))
  output_file_name <- paste(p_file_name, d_file_name, sep = "-", collapse = NULL)
  output_file_path <- file.path(output_dir, paste(output_file_name, "csv", sep = "."))
  
  # Save results
  write.csv(res_df, output_file_path, row.names = FALSE)
}