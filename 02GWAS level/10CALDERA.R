# Load required packages
library(Matrix)
library(glmnet)
library(data.table)

# Set paths
caldera_path <- "./caldera"
result_path <- "./result"
dir.create(result_path, showWarnings = FALSE)

# Load z_caldera.R
source(file.path(caldera_path, "z_caldera.R"))

# Get all .magma.preds files
pops_files <- list.files(path = "~/S-PoPS", pattern = "\\.magma\\.preds$", full.names = TRUE)

# Process each .magma.preds file
for (pops_file in pops_files) {
  file_name <- basename(pops_file)
  base_name <- sub("\\.magma\\.preds$", "", file_name)
  magma_file <- file.path("./MAGMA", paste0(base_name, ".magma.genes.out"))
  cs_file <- file.path("./CSGWAS", paste0(base_name, ".txt"))
  
  # Run caldera with error handling
  casualdata <- tryCatch({
    caldera(pops_file = pops_file, magma_file = magma_file, cs_file = cs_file, caldera_path = caldera_path)
  }, error = function(e) {
    message(paste("Error processing", base_name, ":", e$message))
    return(NULL)
  })
  
  if (is.null(casualdata)) next
  
  # Filter rows where caldera > 0.5
  filtered_data <- casualdata[casualdata$caldera > 0.5, ]
  
  # Save results
  output_file <- file.path(result_path, paste0(base_name, ".txt"))
  fwrite(filtered_data, file = output_file, sep = "\t")
}