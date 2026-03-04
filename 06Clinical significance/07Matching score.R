library(readxl)
library(dplyr)
library(openxlsx)
library(foreach)
library(doParallel)

pathology_folder <- "TraitEnrichment"
pharmacology_folder <- "Drugpathway"

pathology_files <- list.files(pathology_folder, pattern = "\\.xlsx$", full.names = TRUE)
pharmacology_files <- list.files(pharmacology_folder, pattern = "\\.xlsx$", full.names = TRUE)

calculate_pairing_score <- function(pathology_data, pharmacology_data) {
  pathology_tier1 <- pathology_data$ID[1:20]
  pathology_tier2 <- pathology_data$ID[21:50]
  pharmacology_tier1 <- pharmacology_data$ID[1:20]
  pharmacology_tier2 <- pharmacology_data$ID[21:50]
  
  x <- sum(pathology_tier1 %in% pharmacology_tier1)
  y <- sum(pathology_tier1 %in% pharmacology_tier2)
  z <- sum(pharmacology_tier1 %in% pathology_tier2)
  
  raw_score <- x * 1 + y * 0.5 + z * 0.5
  final_score <- raw_score / 20
  
  return(final_score)
}

num_cores <- 50
cl <- makeCluster(num_cores)
registerDoParallel(cl)

results <- foreach(pathology_file = pathology_files, .combine = rbind, .packages = c("readxl", "dplyr", "openxlsx")) %dopar% {
  local_results <- data.frame(
    Pathology_File = character(),
    Pharmacology_File = character(),
    Pairing_Score = numeric(),
    stringsAsFactors = FALSE
  )
  
  pathology_data <- read.xlsx(pathology_file) %>% as_tibble()
  
  for (pharmacology_file in pharmacology_files) {
    pharmacology_data <- read.xlsx(pharmacology_file) %>% as_tibble()
    
    score <- calculate_pairing_score(pathology_data, pharmacology_data)
    
    pathology_name <- tools::file_path_sans_ext(basename(pathology_file))
    pharmacology_name <- tools::file_path_sans_ext(basename(pharmacology_file))
    
    local_results <- local_results %>%
      add_row(
        Pathology_File = pathology_name,
        Pharmacology_File = pharmacology_name,
        Pairing_Score = score
      )
  }
  
  return(local_results)
}

stopCluster(cl)

output_file <- "Pairing_Scores.xlsx"
write.xlsx(results, output_file, row.names = FALSE)

cat("Pairing scores calculated and saved to:", output_file, "\n")