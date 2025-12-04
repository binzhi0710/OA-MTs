library(dplyr)

data <- read.csv("baselinedata.csv")
colnames(data) <- gsub("^X", "", colnames(data))

covariates <- c("21022", "31", "22189", "21001")
analyze_columns <- c("131868", "131870", "131872", "131874", "131876", "EarlyOA", "AllOA")

continuous_vars <- c("30740", "30750", "4079", "4080", "102", "30760", "30780",
                     "30870", "30690", "WHR")

binary_vars <- setdiff(names(data), c("eid", covariates, continuous_vars, analyze_columns))

data <- data %>%
  mutate(across(all_of(continuous_vars), as.numeric))

data <- data %>%
  mutate(across(all_of(c(binary_vars, analyze_columns)), as.factor))

data <- data %>%
  mutate(
    `21022` = as.numeric(`21022`),
    `31` = as.factor(`31`),
    `22189` = as.numeric(`22189`)
  )

output_dir <- "0909BMIadj"
if (!dir.exists(output_dir)) dir.create(output_dir)

for (outcome in analyze_columns) {
  for (exposure in c(continuous_vars, binary_vars)) {

    covariates_quoted <- sapply(covariates, function(x) paste0("`", x, "`"))
    formula <- as.formula(paste0("`", outcome, "` ~ `", exposure, "` + ",
                                 paste(covariates_quoted, collapse = " + ")))

    model <- glm(formula, data = data, family = binomial)

    output_file <- file.path(output_dir, paste0(outcome, "_", exposure, ".txt"))
    sink(output_file)
    cat("Formula: ", deparse(formula), "\n\n")
    print(summary(model))
    sink()
  }
}
