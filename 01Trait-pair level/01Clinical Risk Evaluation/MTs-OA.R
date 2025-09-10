library(dplyr)

# Load data
data <- read.csv("baselinedata.csv")
colnames(data) <- gsub("^X", "", colnames(data))

# Define covariates (excluding 21001 for BMI analysis)
covariates <- c("21022", "31", "22189", "21001")  # age, sex, Townsend deprivation index, BMI

# Define six target variables (outcomes, binary)
analyze_columns <- c("131868", "131870", "131872", "131874", "131876", "EarlyOA", "AllOA")

# Define continuous variables (exposures)
continuous_vars <- c("30740", "30750", "4079", "4080", "102", "30760", "30780", 
                     "30870", "30690", "WHR")

# Define binary variables (exposures), explicitly excluding analyze_columns
binary_vars <- setdiff(names(data), c("eid", covariates, continuous_vars, analyze_columns))

# Data preprocessing
# Ensure continuous_vars are numeric
data <- data %>%
  mutate(across(all_of(continuous_vars), as.numeric))

# Ensure binary_vars and analyze_columns are binary (0/1), convert to factor
data <- data %>%
  mutate(across(all_of(c(binary_vars, analyze_columns)), as.factor))

# Ensure covariates have correct format
data <- data %>%
  mutate(
    `21022` = as.numeric(`21022`),  # age, ensure numeric
    `31` = as.factor(`31`),         # sex, convert to factor
    `22189` = as.numeric(`22189`)   # Townsend index, ensure numeric
  )

# Check if variables exist
cat("\nChecking if variables exist:\n")
cat("Are all continuous_vars in data: ", all(continuous_vars %in% colnames(data)), "\n")
cat("Are all binary_vars in data: ", all(binary_vars %in% colnames(data)), "\n")
cat("Are all analyze_columns in data: ", all(analyze_columns %in% colnames(data)), "\n")
cat("Are all covariates in data: ", all(covariates %in% colnames(data)), "\n")

# Count number of samples with value 1 for each outcome variable
cat("\nBaseline prevalence count (number of samples with value 1):\n")
for (col in analyze_columns) {
  count_ones <- sum(data[[col]] == "1", na.rm = TRUE)
  cat("Column ", col, " number of 1s: ", count_ones, "\n")
}

# Create output folder (if it doesn't exist)
output_dir <- "0909BMIadj"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop through each outcome variable in analyze_columns and each exposure variable
for (outcome in analyze_columns) {
  cat("\n\n=== Model results for outcome variable ", outcome, " ===\n")
  
  # Loop through continuous_vars and binary_vars
  for (exposure in c(continuous_vars, binary_vars)) {
    # Build formula, adding backticks for numeric column names
    covariates_quoted <- sapply(covariates, function(x) paste0("`", x, "`"))
    formula <- as.formula(paste0("`", outcome, "` ~ `", exposure, "` + ", 
                                 paste(covariates_quoted, collapse = " + ")))
    
    # Print formula for debugging
    cat("\nFormula: ", deparse(formula), "\n")
    
    # Fit logistic regression model
    tryCatch({
      model <- glm(formula, data = data, family = binomial)
      
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