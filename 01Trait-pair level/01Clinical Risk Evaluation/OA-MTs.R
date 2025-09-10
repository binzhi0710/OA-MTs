library(dplyr)

# Load data
data <- read.csv("filtered_data4_transformed.csv")
colnames(data) <- gsub("^X", "", colnames(data))

# Define covariates (excluding 21001 for BMI analysis)
covariates <- c("21022", "31", "22189", "21001")  # age, sex, Townsend deprivation index, BMI

# Define new exposure variables (former outcome variables, binary)
new_exposures <- c("131868", "131870", "131872", "131874", "131876", "EarlyOA", "AllOA")

# Define new outcome variables (former exposure variables)
continuous_outcomes <- c("30740", "30750", "4079", "4080", "102", "30760", "30780", 
                         "30870", "30690", "48", "49", "21001")  # continuous outcomes
binary_outcomes <- setdiff(names(data), c("eid", covariates, continuous_outcomes, new_exposures))  # binary outcomes

# Data preprocessing
# Ensure continuous_outcomes are numeric
data <- data %>%
  mutate(across(all_of(continuous_outcomes), as.numeric))

# Ensure binary_outcomes and new_exposures are binary (0/1), convert to factor
data <- data %>%
  mutate(across(all_of(c(binary_outcomes, new_exposures)), as.factor))

# Ensure covariates have correct format
data <- data %>%
  mutate(
    `21022` = as.numeric(`21022`),  # age, ensure numeric
    `31` = as.factor(`31`),         # sex, convert to factor
    `22189` = as.numeric(`22189`),  # Townsend index, ensure numeric
    `21001` = as.numeric(`21001`)   # BMI, ensure numeric
  )

# Check if variables exist
cat("\nChecking if variables exist:\n")
cat("Are all continuous_outcomes in data: ", all(continuous_outcomes %in% colnames(data)), "\n")
cat("Are all binary_outcomes in data: ", all(binary_outcomes %in% colnames(data)), "\n")
cat("Are all new_exposures in data: ", all(new_exposures %in% colnames(data)), "\n")
cat("Are all covariates in data: ", all(covariates %in% colnames(data)), "\n")

# Count number of samples with value 1 for each binary outcome variable (continuous outcomes excluded)
cat("\nBaseline prevalence count (number of samples with value 1, binary outcomes):\n")
for (col in binary_outcomes) {
  count_ones <- sum(data[[col]] == "1", na.rm = TRUE)
  cat("Column ", col, " number of 1: ", count_ones, "\n")
}

# Create output folder (if it doesn't exist)
output_dir <- "0909noBMIreverse"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop through each new outcome variable and each new exposure variable
all_outcomes <- c(continuous_outcomes, binary_outcomes)
for (outcome in all_outcomes) {
  cat("\n\n=== Model results for outcome variable ", outcome, " ===\n")
  
  # Determine model type
  family_type <- if (outcome %in% continuous_outcomes) "gaussian" else "binomial"
  
  # Loop through new_exposures
  for (exposure in new_exposures) {
    # Build formula, adding backticks for numeric column names
    covariates_quoted <- sapply(covariates, function(x) paste0("`", x, "`"))
    formula <- as.formula(paste0("`", outcome, "` ~ `", exposure, "` + ", 
                                 paste(covariates_quoted, collapse = " + ")))
    
    # Print formula for debugging
    cat("\nFormula: ", deparse(formula), "\n")
    
    # Fit regression model
    tryCatch({
      model <- glm(formula, data = data, family = family_type)
      
      # Print model results to console
      cat("\n--- Model for exposure variable ", exposure, " ---\n")
      print(summary(model))
      
      # Save model results to txt file
      output_file <- file.path(output_dir, paste0(outcome, "_", exposure, ".txt"))
      sink(output_file)
      cat("Formula: ", deparse(formula), "\n\n")
      cat("Model results for outcome variable ", outcome, " and exposure variable ", exposure, "\n")
      print(summary(model))
      sink()
      
    }, error = function(e) {
      cat("\nError in model for outcome ", outcome, " and exposure ", exposure, ": ", e$message, "\n")
    })
  }
}