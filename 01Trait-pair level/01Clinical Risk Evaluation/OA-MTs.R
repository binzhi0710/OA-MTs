library(dplyr)

data <- read.csv("filtered_data4_transformed.csv")
colnames(data) <- gsub("^X", "", colnames(data))

covariates <- c("21022", "31", "22189", "21001")
new_exposures <- c("131868", "131870", "131872", "131874", "131876", "EarlyOA", "AllOA")

continuous_outcomes <- c("30740", "30750", "4079", "4080", "102", "30760", "30780", 
                         "30870", "30690", "48", "49", "21001")

binary_outcomes <- setdiff(names(data), c("eid", covariates, continuous_outcomes, new_exposures))

data <- data %>%
  mutate(across(all_of(continuous_outcomes), as.numeric))

data <- data %>%
  mutate(across(all_of(c(binary_outcomes, new_exposures)), as.factor))

data <- data %>%
  mutate(
    `21022` = as.numeric(`21022`),
    `31` = as.factor(`31`),
    `22189` = as.numeric(`22189`),
    `21001` = as.numeric(`21001`)
  )

output_dir <- "0909noBMIreverse"
if (!dir.exists(output_dir)) dir.create(output_dir)

all_outcomes <- c(continuous_outcomes, binary_outcomes)

for (outcome in all_outcomes) {
  
  family_type <- if (outcome %in% continuous_outcomes) "gaussian" else "binomial"
  
  for (exposure in new_exposures) {
    
    covariates_quoted <- sapply(covariates, function(x) paste0("`", x, "`"))
    formula <- as.formula(paste0("`", outcome, "` ~ `", exposure, "` + ",
                                 paste(covariates_quoted, collapse = " + ")))
    
    model <- glm(formula, data = data, family = family_type)
    
    output_file <- file.path(output_dir, paste0(outcome, "_", exposure, ".txt"))
    sink(output_file)
    cat("Formula: ", deparse(formula), "\n\n")
    print(summary(model))
    sink()
  }
}
