library(readxl)
library(glmnet)
library(caret)
library(e1071)
library(pROC)

expression_data <- read.table("trainingsample_z_score_normalized_data.txt", header = TRUE, row.names = 1)
group_data <- read_excel("group_info.xlsx")
group_data$type <- factor(group_data$type, levels = c("control", "OA"))
group_data$name <- as.character(group_data$name)

valid_samples <- intersect(group_data$name, colnames(expression_data))
group_data <- group_data[group_data$name %in% valid_samples, ]
expression_data <- expression_data[, valid_samples, drop = FALSE]

allCgenes <- readLines("allCgenes.txt")
filtered_data <- expression_data[rownames(expression_data) %in% allCgenes, , drop = FALSE]
prepared_data <- t(filtered_data)

train_labels <- ifelse(group_data$type == "control", 0, 1)
input <- cbind(Label = train_labels, as.data.frame(prepared_data))

source("msvmRFE.R")
svmRFE(input, k = 5, halve.above = 100)

nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len = nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k = 5, halve.above = 100)

top.features = WriteFeatures(results, input, save = FALSE)
head(top.features)
write.csv(top.features, "feature_svm_top_features.csv")

ncol = ncol(input)
featseep = lapply(1:ncol, FeatSweep.wrap, results, input)
save(featseep, file = "svm_result.RData")

no.info = min(prop.table(table(input[,1])))
errors = sapply(featseep, function(x) ifelse(is.null(x), NA, x$error))

PlotErrors(errors, no.info = no.info)

which.min(errors)
top <- top.features[1:which.min(errors), "FeatureName"]
write.csv(top, "svm_selected_core_genes.csv")