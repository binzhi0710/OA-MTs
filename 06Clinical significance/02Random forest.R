library(readxl)
library(glmnet)
library(caret)
library(e1071)
library(pROC)
library(randomForest)

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

valid_samples <- intersect(group_data$name, rownames(prepared_data))
group_data <- group_data[group_data$name %in% valid_samples, ]
prepared_data <- prepared_data[valid_samples, , drop = FALSE]

print(dim(prepared_data))
print(dim(group_data))

rf_data <- as.data.frame(prepared_data)
colnames(rf_data) <- make.names(colnames(rf_data))
rf_data$Group <- group_data$type

set.seed(1234)
rf <- randomForest(Group ~ ., data = rf_data, ntree = 500, importance = TRUE)

pdf(file = "RF_Error_Rates.pdf", width = 6, height = 6)
plot(rf, main = "Random Forest Error Rates", lwd = 2)
dev.off()

optimal_trees <- which.min(rf$err.rate[, "OOB"])
cat("Optimal number of trees:", optimal_trees, "\n")

rf_opt <- randomForest(Group ~ ., data = rf_data, ntree = optimal_trees)

importance_rf <- importance(rf_opt)
importance_df <- as.data.frame(importance_rf)
importance_df$Gene <- rownames(importance_rf)
importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]

write.csv(importance_df, "RF_Gene_Importance_Scores.csv", row.names = FALSE)