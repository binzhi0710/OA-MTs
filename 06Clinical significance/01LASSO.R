# Read expression matrix
expression_data <- read.table("trainingsample_z_score_normalized_data.txt", header = TRUE, row.names = 1)

# Read group information
group_data <- read_excel("group_info.xlsx")

# Convert group information types
group_data$type <- factor(group_data$type, levels = c("control", "OA"))
group_data$name <- as.character(group_data$name)  # Ensure sample names are character vectors

# Filter valid samples
valid_samples <- intersect(group_data$name, colnames(expression_data))
group_data <- group_data[group_data$name %in% valid_samples, ]
expression_data <- expression_data[, valid_samples, drop = FALSE]

# 2. Load `allCgenes` gene list
allCgenes <- readLines("allCgenes.txt")

# Filter gene expression data for `allCgenes`
filtered_data <- expression_data[rownames(expression_data) %in% allCgenes, , drop = FALSE]

# Transpose expression matrix so rows are samples and columns are genes
prepared_data <- t(filtered_data)

# Create classification labels vector, control = 0, OA = 1
train_labels <- ifelse(group_data$type == "control", 0, 1)

# Save expression matrix and classification labels as CSV files
write.csv(prepared_data, "expression_matrix.csv", row.names = TRUE)
write.csv(train_labels, "classification_labels.csv", row.names = TRUE)

# 3. Load and convert to matrix format
x <- as.matrix(read.csv("expression_matrix.csv", row.names = 1))  # Expression matrix
y <- as.matrix(read.csv("classification_labels.csv", row.names = 1))  # Classification labels
x <- as.matrix(x)
y <- as.matrix(y)

set.seed(1234)

cvfit <- cv.glmnet(x, y, type.measure = "mse", nfolds = 5, alpha = 1)

# Plot
png("LASSO_CV_Results.png", width = 800, height = 600)
plot(cvfit)
title("LASSO Cross-Validation Results")
dev.off()

# Extract the minimum lambda value
lambda_min <- cvfit$lambda.min
lambda_1se <- cvfit$lambda.1se
cat("Minimum lambda value:", lambda_min, "\n")
cat("1-SE lambda value:", lambda_1se, "\n")

# Rebuild LASSO model using the optimal lambda value
lasso <- glmnet(x, y, family = "binomial", alpha = 1, nlambda = 100)

# Plot LASSO coefficient path and save the image
png("LASSO_Coefficient_Path.png", width = 800, height = 600)
plot(lasso, xvar = "lambda", label = TRUE)
title("LASSO Coefficient Path")
dev.off()

coef <- coef(lasso, s = lambda_min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene <- row.names(coef)[index]

geneCoef <- cbind(Gene = lassoGene, Coef = actCoef)
write.csv(geneCoef, "core_gene_coefficients.csv", row.names = FALSE)
