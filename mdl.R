library(ggplot2)
setwd("C:/Users/JIC/Documents/UNIPD/Stat Learning/proj")
omics <- read.table("data.csv", header = TRUE, sep=",")
rs.cols <- grep(paste0("^", "rs"), names(omics))
cn.cols <- grep(paste0("^", "cn"), names(omics))
mu.cols <- grep(paste0("^", "mu"), names(omics))
pp.cols <-grep(paste0("^", "pp"), names(omics))

# Missing value imputation
# due to high dim, removing missing not an option; every row and col has missing values
rs <- omics[rs.cols]
pp <- omics[pp.cols]
rs.missing <- sum(apply(omics[rs.cols], 1, function(x) any(x == 0))) 
rs.medians <- apply(rs, 2, function(x) median(x[x != 0]))
pp.medians <- apply(pp, 2, function(x) median(x[x != 0]))
# Missing value rate
rs.missing.rate <- sum(rs == 0) / (length(rs.cols) * length(rs))
rs.missing.rate
pp.missing.rate <- sum(pp == 0) / (length(pp.cols) * length(pp))
pp.missing.rate

for(i in seq_along(rs.medians)) {
  rs[rs[, i] == 0, i] <- rs.medians[i]
}
for(i in seq_along(pp.medians)) {
  pp[pp[, i] == 0, i] <- pp.medians[i]
}
omics[c(rs.cols,pp.cols)] <- c(rs,pp)

# Feature distributions
# Gene Expression
long_rs <- tidyr::gather(omics[rs.cols], value = "value", rs.cols)
ggplot(long_rs, aes(x = value)) + 
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Smoothed Distribution of Gene Expression Levels", x = "Level", y = "Density")
# Protein Abundance
long_pp <- tidyr::gather(omics[pp.cols], value = "value", pp.cols)
ggplot(long_pp, aes(x = value)) + 
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Smoothed Distribution of Protein Abundance Levels", x = "Level", y = "Density")
# Copy Number
long_cn <- tidyr::gather(omics[cn.cols], value = "value", cn.cols)$value
barplot(table(long_cn)/sum(table(long_cn)))
# Mutations
long_mu <- tidyr::gather(omics[mu.cols], value = "value", mu.cols)$value
barplot(table(long_mu)/sum(table(long_mu)))
# Vital Status (Target)
barplot(table(omics$vital.status)/sum(table(omics$vital.status)))

# Normalize Continuous Features (Expression and Protein Abundance)
omics_norm <- omics
omics_norm[c(rs.cols,pp.cols)] <- scale(omics[c(rs.cols,pp.cols)])
colMeans(omics)[1:5]

# Ridge Regression
library(glmnet)
X <- omics_norm[,-ncol(omics_norm)]
y <- omics_norm$vital.status 
glm.out <- glm(vital.status~., data=omics_norm, family=binomial)

grid <- 10^seq(10, -2, length=100)
lasso.mod <- glmnet(X,y,alpha=1)
plot(lasso.mod, xvar="lambda", label=TRUE)
plot(lasso.mod, xvar="norm", label=TRUE)
plot(lasso.mod, xvar="dev", label=TRUE)


# TRAIN TEST SPLIT TEST
library(glmnet)

# Assuming omics_norm is your dataset
X <- as.matrix(omics_norm[,-ncol(omics_norm)]) # Convert to matrix for glmnet
y <- omics_norm$vital.status

# Split the data into training (70%) and testing (30%) sets
set.seed(123) # for reproducibility
sample_size <- floor(0.7 * nrow(omics_norm))
train_indices <- sample(seq_len(nrow(omics_norm)), size = sample_size)
X_train <- X[train_indices, ]
y_train <- y[train_indices]
X_test <- X[-train_indices, ]
y_test <- y[-train_indices]

# Fit the model on the training set
grid <- 10^seq(10, -2, length = 100)
lasso.mod <- glmnet(X_train, y_train, alpha = 1, lambda = grid)

# Use cross-validation to find the optimal lambda value
cv.lasso <- cv.glmnet(X_train, y_train, alpha = 1, type.measure = "class")

# Use the best lambda from cross-validation
best.lambda <- cv.lasso$lambda.min

# Make predictions on the test set
predictions <- predict(cv.lasso, newx = X_test, s = "lambda.min", type = "response")
predicted.classes <- ifelse(predictions > 0.5, 1, 0)

# Calculate accuracy
correct_predictions <- sum(predicted.classes == y_test)
accuracy <- correct_predictions / length(y_test)
print(paste("Accuracy on test set:", accuracy))


# Get number non-zero coeffs
# Extract coefficients at the best lambda value
coefficients <- coef(cv.lasso, s = "lambda.min")

# The 'coefficients' object includes a column for the intercept. To count only features:
# Exclude the intercept by starting from the second element
non_zero_features <- coefficients[-1, , drop = FALSE]

# Count the number of non-zero coefficients (features)
num_non_zero_features <- sum(non_zero_features != 0)

print(paste("Number of non-zero features used in the model:", num_non_zero_features))

## Get the names of the features
coef_vector <- as.vector(coefficients[-1, , drop = FALSE])

# Get the names of all features
feature_names <- rownames(coefficients)[-1] # Exclude the intercept row

# Identify the names of the features with non-zero coefficients
non_zero_feature_names <- feature_names[coef_vector != 0]

print("Features with non-zero coefficients:")
print(non_zero_feature_names)
